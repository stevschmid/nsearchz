const std = @import("std");
const print = std.debug.print;
const ArrayList = std.ArrayList;
const Allocator = std.mem.Allocator;

const utils = @import("utils.zig");

const Sequence = @import("sequence.zig").Sequence;
const FastaReader = @import("fasta_reader.zig").FastaReader;
const Database = @import("database.zig").Database;

const bio = @import("bio/bio.zig");
const alphabet = bio.alphabet;

const alphabetChosen = alphabet.DNA;
const kmerInfo = bio.kmer.KmerInfo(alphabetChosen, 8);

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    var bench_start = std.time.milliTimestamp();
    var reader = FastaReader(alphabetChosen).init(allocator);
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
    print("Reading took {}ms\n", .{std.time.milliTimestamp() - bench_start});

    var db = try Database(alphabet.DNA, 8).init(allocator, reader.sequences.items);
    defer db.deinit();
}
