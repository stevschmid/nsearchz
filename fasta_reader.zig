const std = @import("std");
const Sequence = @import("sequence.zig").Sequence;
const SequenceList = @import("sequence.zig").SequenceList;

const utils = @import("utils.zig");

pub fn FastaReader(comptime A: type) type {
    return struct {
        pub fn parse(reader: anytype, sequences: *SequenceList(A)) !void {
            var buffered_reader = std.io.bufferedReader(reader);
            var stream = buffered_reader.reader();

            var identifier = std.ArrayList(u8).init(sequences.allocator);
            defer identifier.deinit();

            var data = std.ArrayList(u8).init(sequences.allocator);
            defer data.deinit();

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
                        // ensure upper case
                        for (line) |*letter| {
                            if (letter.* >= 'a' and letter.* <= 'z') {
                                letter.* -= ('a' - 'A');
                            }
                        }

                        try data.appendSlice(line);
                    },
                }
            }

            // add remaining sequence
            if (data.items.len > 0) {
                try sequences.append(identifier.items, data.items);
            }
        }

        pub fn readFile(path: []const u8, sequences: *SequenceList(A)) !void {
            const dir: std.fs.Dir = std.fs.cwd();
            const file: std.fs.File = try dir.openFile(path, .{ .read = true });
            defer file.close();

            return try parse(file.reader(), sequences);
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

    var sequences = SequenceList(DNA).init(allocator);
    defer sequences.deinit();

    var source = std.io.StreamSource{ .const_buffer = std.io.fixedBufferStream(fasta) };
    try FastaReader(DNA).parse(source.reader(), &sequences);

    try std.testing.expectEqual(@as(usize, 3), sequences.list.items.len);

    try std.testing.expectEqualStrings("Seq1", sequences.list.items[0].identifier);
    try std.testing.expectEqualStrings("TGGCGAAATTGGG", sequences.list.items[0].data);

    try std.testing.expectEqualStrings("Seq2", sequences.list.items[1].identifier);
    try std.testing.expectEqualStrings("TTTTTCAGTC", sequences.list.items[1].data);

    try std.testing.expectEqualStrings("Seq3", sequences.list.items[2].identifier);
    try std.testing.expectEqualStrings("ACTGC", sequences.list.items[2].data);
}
