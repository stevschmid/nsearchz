const std = @import("std");
const print = std.debug.print;
const ArrayList = std.ArrayList;
const Allocator = std.mem.Allocator;

const DNA = struct {
    pub fn complement(letter: u8) u8 {
        _ = letter;
        return 'Y';
    }
};

const Protein = struct {
    pub fn complement(letter: u8) u8 {
        _ = letter;
        return 'X';
    }
};

pub fn dup(comptime T: type, allocator: Allocator, slice: []const T) ![]T {
    const result = try allocator.alloc(T, slice.len);
    std.mem.copy(u8, result, slice);
    return result;
}

pub fn Sequence(comptime A: type) type {
    return struct {
        const Self = @This();

        allocator: Allocator,

        identifier: []const u8,
        data: []const u8,

        pub fn init(allocator: Allocator, identifier: []const u8, data: []const u8) !Self {
            return Self {
                .allocator = allocator,
                .identifier = try dup(u8, allocator, identifier),
                .data = try dup(u8, allocator, data),
            };
        }

        pub fn deinit(self: *Self) void {
            self.allocator.free(self.identifier);
            self.allocator.free(self.data);
        }

        pub fn matches(self: *Self, other: Self) bool {
            if (self.data.len != other.data.len) {
                return false;
            }

            return for (self.data) |_, index| {
                if (self.data[index] != other.data[index]) break false;
            } else true;
        }

        pub fn complement(self: *Self) !Self {
            var compl_data = try dup(u8, self.allocator, self.data);

            for (compl_data) |*letter| {
                letter.* = A.complement(letter.*);
            }

            return Sequence(A) {
                .allocator = self.allocator,
                .identifier = try dup(u8, self.allocator, self.identifier),
                .data = compl_data,
            };
        }
    };
}

// pub fn addSequence(list:*std.ArrayList(Sequence), allocator: Allocator, identifier:[]const u8, data:[]const u8) !void {
//     if (data.len <= 0) return;

//     try list.append(Sequence{ 
//         .allocator = allocator,
//         .identifier = identifier,
//         .data = data,
//     });
// }

pub fn FastaReader(comptime A: type) type {
    return struct {
        const Self = @This();

        allocator: std.mem.Allocator,
        sequences: std.ArrayList(Sequence(A)),

        pub fn init(allocator: Allocator) Self {
            return Self{
                .allocator = allocator,
                .sequences = std.ArrayList(Sequence(A)).init(allocator),
            };
        }

        pub fn deinit(self: *Self) void {
            for (self.sequences.items) |*sequence| sequence.deinit();
            self.sequences.deinit();
        }

        pub fn readFile(self: *Self, path: []const u8) !void {
            _ = self;
            _ = path;

            const dir: std.fs.Dir = std.fs.cwd();
            const file: std.fs.File = try dir.openFile(path, .{ .read = true });
            defer file.close();

            var buffered_reader = std.io.bufferedReader(file.reader());
            var stream = buffered_reader.reader();

            var identifier = std.ArrayList(u8).init(self.allocator);
            defer identifier.deinit();

            var data = std.ArrayList(u8).init(self.allocator);
            defer data.deinit();

            while (try stream.readUntilDelimiterOrEofAlloc(self.allocator, '\n', 1024)) |line| {
                defer self.allocator.free(line);

                switch (line[0]) {
                    '>' => {
                        print("New {s}\n", .{line});
                        try self.addSequence(identifier.items, data.items);

                        // reset
                        identifier.clearRetainingCapacity();
                        data.clearRetainingCapacity();

                        // start new sequence
                        try identifier.appendSlice(line[1..]);
                    },
                    ';' => {
                        print("Comment {s}\n", .{line});
                        // ignore
                    },
                    else => {
                        print("Data {s} {}\n", .{line, line.len});
                        // data
                        try data.appendSlice(line);
                    }
                }
            }
            // try self.addSequence(identifier.items, data.items);
        }

        fn addSequence(self: *Self, identifier: []const u8, data: []const u8) !void {
            if (data.len == 0) {
                return;
            }

            var sequence = try Sequence(A).init(self.allocator, identifier, data);
            try self.sequences.append(sequence);
        }
    };
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    var reader = FastaReader(DNA).init(allocator);
    defer reader.deinit();

    try reader.readFile("test.fasta");

    // // var compl = try seq.complement();
    // // defer compl.deinit();

    // // print("Compl %{s}\n", .{compl.data});

    // // var comp: Sequence(DNA) = seq.

    // // var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    // // defer _ = gpa.deinit();
    // // const allocator = gpa.allocator();

    // const dir: std.fs.Dir = std.fs.cwd();
    // const file: std.fs.File = try dir.openFile("test.fasta", .{ .read = true });
    // defer file.close();

    // var buffered_reader = std.io.bufferedReader(file.reader());
    // var stream = buffered_reader.reader();

    // var sequences = ArrayList(Sequence(DNA)).init(allocator);
    // defer {
    //     for (sequences.items) |*seq| seq.deinit();
    //     sequences.deinit();
    // }

    // var identifier = ArrayList(u8).init(allocator);
    // defer identifier.deinit();

    // var data = ArrayList(u8).init(allocator);
    // defer data.deinit();

    // // var comment = ArrayList(u8).init(allocator);
    // // defer comment.deinit();

    // while (try stream.readUntilDelimiterOrEofAlloc(allocator, '\n', 1024)) |line| {
    //     defer allocator.free(line);

    //     switch (line[0]) {
    //         '>' => {
    //             // print("New {s}\n", .{line});
    //             try addSequence(&sequences, allocator, identifier.items, data.items);
    //             identifier.clearAndFree();
    //             data.clearAndFree();
    //             try identifier.appendSlice(line[1..]);
    //         },
    //         ';' => {
    //             // print("Comment {s}\n", .{line});
    //             // try comment.appendSlice(line[1..]);
    //         },
    //         else => {
    //             // print("Data {s} {}\n", .{line, line.len});
    //             try data.appendSlice(line);
    //         }
    //     }
    // }

    // // for (sequences.items) |seq| {
    // //     print("Sequence {s}\n", .{seq.identifier});
    // //     print("Data\n", .{});
    // //     print("{s}\n", .{seq.data});
    // //     print("\n", .{});
    // // }
}
