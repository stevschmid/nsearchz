const std = @import("std");
const print = std.debug.print;
const ArrayList = std.ArrayList;
const Allocator = std.mem.Allocator;

const DNA = @import("alphabet/dna.zig").DNA;
const Protein = @import("alphabet/protein.zig").Protein;

const dup = @import("utils.zig").dup;

const Sequence = @import("sequence.zig").Sequence;
const FastaReader = @import("fasta_reader.zig").FastaReader;

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
        std.debug.print("Expected first argument to be path to input file\n", .{});
        return error.InvalidArgs;
    });
    defer allocator.free(file);

    try reader.readFile(file);

    for (reader.sequences.items) |sequence| {
        print("{s}\n{s}\n\n", .{sequence.identifier, sequence.data});

        var compl = try sequence.complementAlloc();
        defer compl.deinit();

        print("Complement\n{s}\n{s}\n\n", .{compl.identifier, compl.data});
    }
}
