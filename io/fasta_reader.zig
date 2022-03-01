const std = @import("std");

const Sequence = @import("../sequence.zig").Sequence;
const utils = @import("../utils.zig");

pub fn FastaReader(comptime A: type) type {
    return struct {
        const Self = @This();

        allocator: std.mem.Allocator,

        stream: std.io.PeekStream(.{ .Static = 4096 }, std.io.BufferedReader(4096, std.io.StreamSource.Reader)),

        identifier: std.ArrayList(u8),
        data: std.ArrayList(u8),

        pub fn init(allocator: std.mem.Allocator, source: *std.io.StreamSource) Self {
            return Self{
                .allocator = allocator,

                .stream = std.io.peekStream(4096, std.io.bufferedReader(source.reader())),

                .identifier = std.ArrayList(u8).init(allocator),
                .data = std.ArrayList(u8).init(allocator),
            };
        }

        pub fn deinit(self: *Self) void {
            self.identifier.deinit();
            self.data.deinit();
        }

        pub fn next(self: *Self) !?Sequence(A) {
            var buffer: [4096]u8 = undefined;

            var reader = self.stream.reader();

            while (try reader.readUntilDelimiterOrEof(&buffer, '\n')) |line| {
                if (line.len == 0)
                    continue;

                switch (line[0]) {
                    '>' => {
                        if (self.identifier.items.len > 0) {
                            // seek back, this line is start of next sequence
                            try self.stream.putBackByte('\n');
                            try self.stream.putBack(line);

                            return self.publish();
                        } else {
                            try self.identifier.appendSlice(line[1..]);
                        }
                    },
                    ';' => {
                        // ignore comment
                    },
                    else => {
                        // ensure upper case
                        for (line) |*letter| {
                            if (letter.* >= 'a' and letter.* <= 'z') {
                                letter.* -= ('a' - 'A');
                            }
                        }

                        try self.data.appendSlice(line);
                    },
                }
            }

            return self.publish();
        }

        fn publish(self: *Self) !?Sequence(A) {
            if (self.identifier.items.len == 0)
                return null;

            if (self.data.items.len == 0)
                return null;

            // build this sequence
            const sequence = try Sequence(A).init(self.allocator, self.identifier.items, self.data.items);

            // reset
            self.identifier.clearRetainingCapacity();
            self.data.clearRetainingCapacity();

            return sequence;
        }
    };
}

const DNA = @import("../bio/bio.zig").alphabet.DNA;

test "reads fasta" {
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

    var reader = FastaReader(DNA).init(allocator, &source);
    defer reader.deinit();

    var seq: Sequence(DNA) = undefined;

    seq = (try reader.next()).?;
    try std.testing.expectEqualStrings(seq.identifier, "Seq1");
    try std.testing.expectEqualStrings(seq.data, "TGGCGAAATTGGG");
    seq.deinit();

    seq = (try reader.next()).?;
    try std.testing.expectEqualStrings(seq.identifier, "Seq2");
    try std.testing.expectEqualStrings(seq.data, "TTTTTCAGTC");
    seq.deinit();

    seq = (try reader.next()).?;
    try std.testing.expectEqualStrings(seq.identifier, "Seq3");
    try std.testing.expectEqualStrings(seq.data, "ACTGC");
    seq.deinit();

    try std.testing.expect((try reader.next()) == null);
    try std.testing.expect((try reader.next()) == null);
}
