const std = @import("std");
const print = std.debug.print;
const ArrayList = std.ArrayList;
const Allocator = std.mem.Allocator;

const utils = @import("utils.zig");
const ea = @import("extend_align.zig");
const ba = @import("banded_align.zig");

const Cigar = @import("cigar.zig").Cigar;
const CigarOp = @import("cigar.zig").CigarOp;

const Sequence = @import("sequence.zig").Sequence;
const SequenceList = @import("sequence.zig").SequenceList;

const FastaReader = @import("fasta_reader.zig").FastaReader;
const Database = @import("database.zig").Database;
const Highscores = @import("highscores.zig").Highscores;

const bio = @import("bio/bio.zig");
const alphabet = bio.alphabet;

const alphabetChosen = alphabet.DNA;
const kmerInfo = bio.kmer.KmerInfo(alphabetChosen, 8);

const Strand = enum {
    Plus,
    Minus,
    Both,
};

const SearchParams = struct {
    max_accepts: u32 = 1,
    max_rejects: u32 = 16,
    min_identity: f32 = 0.8,
    strand: Strand = Strand.Plus,
};

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    var bench_start = std.time.milliTimestamp();

    var arg_it = std.process.args();

    // skip my own exe name
    _ = arg_it.skip();

    const file = (try arg_it.next(allocator) orelse {
        print("Expected first argument to be path to input file\n", .{});
        return error.InvalidArgs;
    });
    defer allocator.free(file);

    var seq_list = SequenceList(alphabetChosen).init(allocator);
    defer seq_list.deinit();

    try FastaReader(alphabetChosen).readFile(file, &seq_list);

    print("Reading took {}ms\n", .{std.time.milliTimestamp() - bench_start});

    bench_start = std.time.milliTimestamp();
    var db = try Database(alphabet.DNA, 8).init(allocator, seq_list.toOwnedSlice());
    defer db.deinit();

    print("Indexing took {}ms\n", .{std.time.milliTimestamp() - bench_start});
}
