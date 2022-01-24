const std = @import("std");
const print = std.debug.print;
const ArrayList = std.ArrayList;
const Allocator = std.mem.Allocator;

const DNA = @import("alphabet/dna.zig").DNA;
const Protein = @import("alphabet/protein.zig").Protein;

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
            var sanitizedData = try dup(u8, allocator, data);

            // convert lower case to upper case
            for (sanitizedData) |*letter| {
                if (letter.* >= 'a' and letter.* <= 'z') {
                    letter.* -= ('a' - 'A');
                }
            }


            return Self {
                .allocator = allocator,
                .identifier = try dup(u8, allocator, identifier),
                .data = sanitizedData,
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

        pub fn read(self: *Self, reader: anytype) !void {
            var buffered_reader = std.io.bufferedReader(reader);
            var stream = buffered_reader.reader();

            var identifier = std.ArrayList(u8).init(self.allocator);
            defer identifier.deinit();

            var data = std.ArrayList(u8).init(self.allocator);
            defer data.deinit();

            var buffer = try self.allocator.alloc(u8, 1024);
            defer self.allocator.free(buffer);

            while (try stream.readUntilDelimiterOrEof(buffer, '\n')) |line| {
                if (line.len == 0) continue;

                switch (line[0]) {
                    '>' => {
                        // add old sequence
                        try self.add(identifier.items, data.items);

                        // reset
                        identifier.clearRetainingCapacity();
                        data.clearRetainingCapacity();

                        // start new sequence
                        try identifier.appendSlice(line[1..]);
                    },
                    ';' => {
                        // ignore comment
                    },
                    else => {
                        try data.appendSlice(line);
                    }
                }
            }
            // add remaining sequence
            try self.add(identifier.items, data.items);
        }

        pub fn readFile(self: *Self, path: []const u8) !void {
            const dir: std.fs.Dir = std.fs.cwd();
            const file: std.fs.File = try dir.openFile(path, .{ .read = true });
            defer file.close();

            try self.read(file.reader());
        }

        fn add(self: *Self, identifier: []const u8, data: []const u8) !void {
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
    }
}

test "reads fasta" {
    const allocator = std.testing.allocator;
    var reader = FastaReader(DNA).init(allocator);
    defer reader.deinit();

    const fasta = 
        \\>Seq1
        \\;comment2
        \\;comment2
        \\TGGCGAA
        \\ATTGGG
        \\
        \\>Seq2
        \\TTTTT
        \\CAGTC
        \\>Seq3
        \\actgc
        ;

    var source = std.io.StreamSource{ .const_buffer = std.io.fixedBufferStream(fasta) };
    try reader.read(source.reader());

    try std.testing.expect(reader.sequences.items.len == 3);

    try std.testing.expectEqualStrings("Seq1", reader.sequences.items[0].identifier);
    try std.testing.expectEqualStrings("TGGCGAAATTGGG", reader.sequences.items[0].data);

    try std.testing.expectEqualStrings("Seq2", reader.sequences.items[1].identifier);
    try std.testing.expectEqualStrings("TTTTTCAGTC", reader.sequences.items[1].data);

    try std.testing.expectEqualStrings("Seq3", reader.sequences.items[2].identifier);
    try std.testing.expectEqualStrings("ACTGC", reader.sequences.items[2].data);
}
