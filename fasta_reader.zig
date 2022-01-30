const std = @import("std");
const Sequence = @import("sequence.zig").Sequence;

pub fn SequenceList(comptime A: type) type {
    return struct {
        const Self = @This();

        allocator: std.mem.Allocator,
        list: std.ArrayList(Sequence(A)),

        pub fn init(allocator: std.mem.Allocator) Self {
            return Self{
                .allocator = allocator,
                .list = std.ArrayList(Sequence(A)).init(allocator),
            };
        }

        pub fn append(self: *Self, identifier: []const u8, data: []const u8) !void {
            var sequence = try Sequence(A).init(self.allocator, identifier, data);
            errdefer sequence.deinit();

            try self.list.append(sequence);
        }

        pub fn deinit(self: *Self) void {
            for (self.list.items) |*seq|
                seq.deinit();

            self.list.deinit();
        }
    };
}

pub fn FastaReader(comptime A: type) type {
    return struct {
        pub fn parseAlloc(allocator: std.mem.Allocator, reader: anytype) !SequenceList(A) {
            var buffered_reader = std.io.bufferedReader(reader);
            var stream = buffered_reader.reader();

            var identifier = std.ArrayList(u8).init(allocator);
            defer identifier.deinit();

            var data = std.ArrayList(u8).init(allocator);
            defer data.deinit();

            var sequences = SequenceList(A).init(allocator);
            errdefer sequences.deinit();

            var buffer: [4096]u8 = undefined;
            while (try stream.readUntilDelimiterOrEof(&buffer, '\n')) |line| {
                if (line.len == 0) continue;

                switch (line[0]) {
                    '>' => {
                        // add old sequence
                        if (data.items.len > 0) {
                            try sequences.append(identifier.items, data.items);
                        }

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
                    },
                }
            }

            // add remaining sequence
            if (data.items.len > 0) {
                try sequences.append(identifier.items, data.items);
            }

            return sequences;
        }

        pub fn parseFileAlloc(allocator: std.mem.Allocator, path: []const u8) ![]Sequence(A) {
            const dir: std.fs.Dir = std.fs.cwd();
            const file: std.fs.File = try dir.openFile(path, .{ .read = true });
            defer file.close();

            return try parseAlloc(allocator, file.reader());
        }
    };
}

test "reads fasta" {
    const DNA = @import("bio/alphabet/dna.zig").DNA;

    const allocator = std.testing.allocator;

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
    var sequences = try FastaReader(DNA).parseAlloc(allocator, source.reader());
    defer sequences.deinit();

    try std.testing.expectEqual(@as(usize, 3), sequences.list.items.len);

    // try std.testing.expectEqualStrings("Seq1", sequences[0].identifier);
    // try std.testing.expectEqualStrings("TGGCGAAATTGGG", sequences[0].data);

    // try std.testing.expectEqualStrings("Seq2", sequences[1].identifier);
    // try std.testing.expectEqualStrings("TTTTTCAGTC", sequences[1].data);

    // try std.testing.expectEqualStrings("Seq3", sequences[2].identifier);
    // try std.testing.expectEqualStrings("ACTGC", sequences[2].data);
}
