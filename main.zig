const std = @import("std");
const print = std.debug.print;
const ArrayList = std.ArrayList;
const Allocator = std.mem.Allocator;

const utils = @import("utils.zig");

const Sequence = @import("sequence.zig").Sequence;
const FastaReader = @import("fasta_reader.zig").FastaReader;

const alphabet = @import("alphabet.zig");

pub fn KmerGenerator(comptime A: type) type {
    return struct {
        const Self = @This();

        const AlphabetInfo = alphabet.AlphabetInfo(A);
        const Mask: AlphabetInfo.KmerType = ~@intCast(AlphabetInfo.KmerType, 0);

        pos: usize,
        letters: []const u8,
        val: AlphabetInfo.KmerType = 0,
        ambiguity_count: usize = 0,

        fn init(letters: []const u8) Self {
            return Self{
                .pos = 0,
                .letters = letters,
            };
        }

        pub fn advance(self: *Self) bool {
            while (self.pos + 1 < AlphabetInfo.NumLettersPerKmer and self.consumeNext()) {
                // advance up to pos - 1 for the initial kmer
            }

            return self.consumeNext();
        }

        pub fn kmer(self: *Self) ?AlphabetInfo.KmerType {
            return if (self.ambiguity_count > 0) null else self.val;
        }

        fn consumeNext(self: *Self) bool {
            if (self.pos >= self.letters.len) {
                return false;
            }

            // evaluate current letter
            const letterBits = A.mapToBits(self.letters[self.pos]);
            if (letterBits == null) {
                // the next X kmers are considered to be ambiguous
                self.ambiguity_count = AlphabetInfo.NumLettersPerKmer + 1;
            } else {
                // map current letter
                self.val = ((self.val << @bitSizeOf(AlphabetInfo.LetterToBitMapType)) | letterBits.?) & AlphabetInfo.KmerMask;
            }

            // advance
            self.pos += 1;

            // in ambiguous region?
            if (self.ambiguity_count > 0)  {
                self.ambiguity_count -= 1;
            }

            return true;
        }
    };
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    // var reader = Reader(Fasta, DNA).init(allocator);

    var reader = FastaReader(alphabet.DNA).init(allocator);
    defer reader.deinit();

    var arg_it = std.process.args();

    // skip my own exe name
    _ = arg_it.skip();

    const file = (try arg_it.next(allocator) orelse {
        print("Expected first argument to be path to input file\n", .{});
        return error.InvalidArgs;
    });
    defer allocator.free(file);
    try reader.readFile(file);

    // build database PogU
    var counts = try allocator.alloc(usize, alphabet.AlphabetInfo(alphabet.DNA).MaxKmers);
    std.mem.set(usize, counts, 0);
    defer allocator.free(counts);

    // indices, for first loop simply to keep track of the unique kmer of a given sequence
    var unique_tracking = try allocator.alloc(isize, alphabet.AlphabetInfo(alphabet.DNA).MaxKmers);
    std.mem.set(isize, unique_tracking, -1);
    defer allocator.free(unique_tracking);

    print("Length {}\n", .{counts.len});

    var total_entries: usize = 0;
    var total_unique_entries: usize = 0;

    for (reader.sequences.items) |sequence, sequence_idx| {
        var kmer_gen = KmerGenerator(alphabet.DNA).init(sequence.data);
        while (kmer_gen.advance()) {
            total_entries += 1;

            const kmer = kmer_gen.kmer();

            // ambiguous?
            if (kmer == null)
                continue;

            // already counted for this sequence?
            if (unique_tracking[kmer.?] == sequence_idx) 
                continue;

            unique_tracking[kmer.?] = @intCast(isize, sequence_idx);
            counts[kmer.?] += 1;
            total_unique_entries += 1;
        }
    }

    var index_offsets = try allocator.alloc(usize, alphabet.AlphabetInfo(alphabet.DNA).MaxKmers);
    defer allocator.free(index_offsets);
    for (index_offsets) |*index_offset, idx| {
        index_offset.* = if (idx > 0) index_offsets[idx - 1] + counts[ idx - 1 ] else 0;
        print("Offset {}\n", .{index_offset.*});
    }

    for (counts) |count, kmer| {
        if (count > 0) {
            print("Kmer {b:0>16}: {}\n", .{kmer, count});
        }
    }
}
