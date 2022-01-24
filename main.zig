const std = @import("std");
const print = std.debug.print;
const ArrayList = std.ArrayList;
const Allocator = std.mem.Allocator;

const Sequence = struct {
    allocator: Allocator,
    sequence: []const u8,
    identifier: []const u8,

    pub fn deinit(self: Sequence) void { 
        self.allocator.free(self.identifier);
        self.allocator.free(self.sequence);
    }
};

pub fn addSequence(list:*std.ArrayList(Sequence), allocator: Allocator, identifier:[]const u8, sequence:[]const u8) !void {
    if (sequence.len <= 0) return;

    try list.append(Sequence{ 
        .allocator = allocator,
        .identifier = identifier,
        .sequence = sequence,
    });
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const dir: std.fs.Dir = std.fs.cwd();
    const file: std.fs.File = try dir.openFile("test.fasta", .{ .read = true });
    defer file.close();

    var buffered_reader = std.io.bufferedReader(file.reader());
    var stream = buffered_reader.reader();

    var sequences = ArrayList(Sequence).init(allocator);
    defer {
        for (sequences.items) |seq| seq.deinit();
        sequences.deinit();
    }

    var identifier = ArrayList(u8).init(allocator);
    defer identifier.deinit();

    var sequence = ArrayList(u8).init(allocator);
    defer sequence.deinit();

    var comment = ArrayList(u8).init(allocator);
    defer comment.deinit();

    while (try stream.readUntilDelimiterOrEofAlloc(allocator, '\n', 1024)) |line| {
        defer allocator.free(line);

        switch (line[0]) {
            '>' => {
                // print("New {s}\n", .{line});
                try addSequence(&sequences, allocator, identifier.toOwnedSlice(), sequence.toOwnedSlice());
                try identifier.appendSlice(line[1..]);
            },
            ';' => {
                // print("Comment {s}\n", .{line});
                try comment.appendSlice(line[1..]);
            },
            else => {
                // print("Data {s} {}\n", .{line, line.len});
                try sequence.appendSlice(line);
            }
        }
    }

    for (sequences.items) |seq| {
        print("Sequence {s}\n", .{seq.identifier});
        print("Data\n", .{});
        print("{s}\n", .{seq.sequence});
        print("\n", .{});
    }
}
