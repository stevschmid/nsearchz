const std = @import("std");
const print = std.debug.print;
const ArrayList = std.ArrayList;
const Allocator = std.mem.Allocator;

const DNA = @import("alphabet/dna.zig").DNA;
const Protein = @import("alphabet/protein.zig").Protein;

const dup = @import("utils.zig").dup;

const Sequence = @import("sequence.zig").Sequence;
const FastaReader = @import("fasta_reader.zig").FastaReader;

pub fn KmerGenerator(comptime A: type, comptime NumLettersPerKmer: comptime_int) type {
    return struct {
        const Self = @This();

        const LetterToBitMapType = @typeInfo(@typeInfo(@TypeOf(DNA.mapToBits)).Fn.return_type.?).Optional.child;
        const NumBitsPerLetter = @bitSizeOf(LetterToBitMapType);

        const NumBitsPerKmer = NumLettersPerKmer * NumBitsPerLetter;

        pub const KmerType = @Type(.{ .Int = .{
            .signedness = .unsigned,
            .bits = NumBitsPerKmer,
        } });

        const Mask: KmerType = ~@intCast(KmerType, 0);

        pos: usize,
        letters: []const u8,
        val: KmerType = 0,
        ambiguity_count: usize = 0,

        fn init(letters: []const u8) Self {
            return Self{
                .pos = 0,
                .letters = letters,
            };
        }

        pub fn advance(self: *Self) bool {
            while (self.pos + 1 < NumLettersPerKmer and self.consumeNext()) {
                // advance up to pos - 1 for the initial kmer
            }

            return self.consumeNext();
        }

        pub fn kmer(self: *Self) ?KmerType {
            return if (self.ambiguity_count > 0) null else self.val;
        }

        fn consumeNext(self: *Self) bool {
            if (self.pos >= self.letters.len) {
                return false;
            }

            // evaluate current letter
            const letterBits = A.mapToBits(self.letters[self.pos]);
            if (letterBits == null) {
                // the next X kmers are considered to be ambigious
                self.ambiguity_count = NumLettersPerKmer + 1;
            } else {
                // map current letter
                self.val = ((self.val << @bitSizeOf(LetterToBitMapType)) | letterBits.?) & Mask;
            }

            // advance
            self.pos += 1;

            // in ambigious region?
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

    var reader = FastaReader(DNA).init(allocator);
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

    for (reader.sequences.items) |sequence| {
        var kmer_gen = KmerGenerator(DNA, 8).init(sequence.data);
        while (kmer_gen.advance()) {
            if (kmer_gen.kmer() != null) {
                print("{b:0>16}\n", .{kmer_gen.kmer()});
            } else  {
                print("Kmer is ambigious\n", .{});
            }
        }
    }
}
