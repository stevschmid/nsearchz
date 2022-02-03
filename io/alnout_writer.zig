const std = @import("std");

const alphabet = @import("../bio/bio.zig").alphabet;
const Sequence = @import("../sequence.zig").Sequence;

const SearchHit = @import("../search.zig").SearchHit;
const SearchHitList = @import("../search.zig").SearchHitList;

const utils = @import("../utils.zig");
const Cigar = @import("../cigar.zig").Cigar;

// Query >QryId
//  %Id   TLen  Target
// 100%     25  DbId3
//  81%     33  DbId1

//  Query 27nt >QryId
// Target 25nt >DbId3

// Qry  1 + GGTGAGACGTTACGCAATAAATTGA 25
//          |||||||||||||||||||||||||
// Tgt  1 + GGTGAGACGTTACGCAATAAATTGA 25

// 25 cols, 25 ids (100.0%), 0 gaps (0.0%)

//  Query 27nt >QryId
// Target 33nt >DbId1

// Qry  1 + GGTGAGACGTTACGCAATAAATTGAGA 27
//           ||||||||  | |||| |||||||||
// Tgt  3 + CGTGAGACG--ATGCAAAAAATTGAGA 27

// 27 cols, 22 ids (81.5%), 2 gaps (7.4%)

pub fn AlnoutWriter(comptime A: type) type {
    return struct {
        pub fn write(writer: anytype, query: Sequence(A), hits: SearchHitList(A)) !void {
            try std.fmt.format(writer, "Query >{s}\n", .{query.identifier});
            try std.fmt.format(writer, " %Id   TLen  Target\n", .{});

            for (hits.list.items) |hit| {
                try std.fmt.format(writer, "{d:3.0}% {d:6}  {s}\n", .{ hit.cigar.identity() * 100.0, hit.target.length(), hit.target.identifier });
            }

            try std.fmt.format(writer, "\n", .{});

            for (hits.list.items) |hit| {
                if (hit.cigar.isEmpty())
                    continue;

                try std.fmt.format(writer, " Query {d}{s} >{s}\n", .{ query.length(), A.Unit, query.identifier });
                try std.fmt.format(writer, "Target {d}{s} >{s}\n", .{ hit.target.length(), A.Unit, hit.target.identifier });

                try std.fmt.format(writer, "\n", .{});

                var query_idx: usize = 0;
                var target_idx: usize = 0;

                var ops = hit.cigar.entries.items;

                const first = ops[0];
                switch (first.op) {
                    .deletion => {
                        target_idx = first.count;
                        ops = ops[1..];
                    },
                    .insertion => {
                        query_idx = first.count;
                        ops = ops[1..];
                    },
                    else => {},
                }

                const last = ops[ops.len - 1];
                switch (last.op) {
                    .deletion => {
                        ops = ops[0 .. ops.len - 1];
                    },
                    .insertion => {
                        ops = ops[0 .. ops.len - 1];
                    },
                    else => {},
                }

                const max_length = std.math.max(query.length(), hit.target.length());

                // Qry  1 + GGTGAGACGTTACGCAATAAATTGA 25
                try std.fmt.format(writer, "Qry ", .{});
                try printLength(writer, max_length, target_idx + 1);

                try std.fmt.format(writer, " + ", .{});

                for (ops) |op| {
                    var count: usize = 0;
                    while (count < op.count) : (count += 1) {
                        var ch = query.data[query_idx];
                        switch (op.op) {
                            .match, .mismatch, .insertion => query_idx += 1,
                            .deletion => ch = '-',
                        }

                        try std.fmt.format(writer, "{c}", .{ch});
                    }
                }

                try std.fmt.format(writer, " ", .{});
                try printLength(writer, max_length, query_idx);
                try std.fmt.format(writer, "\n", .{});

                //          |||||||||||||||||||||||||
                try std.fmt.format(writer, "    ", .{});
                try printPadding(writer, max_length);
                try std.fmt.format(writer, "   ", .{});

                for (ops) |op| {
                    const ch: u8 = switch (op.op) {
                        .match => '|',
                        else => ' ',
                    };

                    var count: usize = 0;
                    while (count < op.count) : (count += 1) {
                        try std.fmt.format(writer, "{c}", .{ch});
                    }
                }

                try std.fmt.format(writer, " ", .{});
                try printPadding(writer, max_length);
                try std.fmt.format(writer, "\n", .{});

                // Tgt  1 + GGTGAGACGTTACGCAATAAATTGA 25
                try std.fmt.format(writer, "Tgt ", .{});
                try printLength(writer, max_length, target_idx + 1);

                try std.fmt.format(writer, " + ", .{});

                for (ops) |op| {
                    var count: usize = 0;
                    while (count < op.count) : (count += 1) {
                        var ch = hit.target.data[target_idx];
                        switch (op.op) {
                            .match, .mismatch, .deletion => target_idx += 1,
                            .insertion => ch = '-',
                        }

                        try std.fmt.format(writer, "{c}", .{ch});
                    }
                }

                try std.fmt.format(writer, " ", .{});
                try printLength(writer, max_length, target_idx);
                try std.fmt.format(writer, "\n", .{});

                try std.fmt.format(writer, "\n", .{});

                // 25 cols, 25 ids (100.0%), 0 gaps (0.0%)
            }
        }

        fn printPadding(writer: anytype, length: usize) !void {
            const max = std.fmt.count("{d}", .{length});
            var cur: usize = 0;

            while (cur < max) : (cur += 1) {
                try std.fmt.format(writer, " ", .{});
            }
        }

        fn printLength(writer: anytype, max_length: usize, length: usize) !void {
            const max = std.fmt.count("{d}", .{max_length});
            var cur = std.fmt.count("{d}", .{length});

            while (cur < max) : (cur += 1) {
                try std.fmt.format(writer, " ", .{});
            }

            try std.fmt.format(writer, "{d}", .{length});
        }
    };
}

test "writes correctly" {
    const DNA = alphabet.DNA;

    const allocator = std.testing.allocator;
    _ = allocator;

    var query = try Sequence(DNA).init(allocator, "QryId", "GGTGAGACGTTACGCAATAAATTGAGA");
    defer query.deinit();

    var cigar = Cigar.init(allocator);
    defer cigar.deinit();

    var hits = SearchHitList(DNA).init(allocator);
    defer hits.deinit();

    // 25=2I
    var target3 = try Sequence(DNA).init(allocator, "DbId3", "GGTGAGACGTTACGCAATAAATTGA");
    defer target3.deinit();

    cigar.clear();
    try cigar.addWithCount(.match, 25);
    try cigar.addWithCount(.insertion, 2);

    try hits.list.append(try SearchHit(DNA).init(allocator, target3, cigar));

    // 2D1X8=2I1=1X4=1X9=1I
    var target1 = try Sequence(DNA).init(allocator, "DbId1", "ATCGTGAGACGATGCAAAAAATTGAGACGGATT");
    defer target1.deinit();

    cigar.clear();
    try cigar.addWithCount(.deletion, 2);
    try cigar.addWithCount(.mismatch, 1);
    try cigar.addWithCount(.match, 8);
    try cigar.addWithCount(.insertion, 2);
    try cigar.addWithCount(.match, 1);
    try cigar.addWithCount(.mismatch, 1);
    try cigar.addWithCount(.match, 4);
    try cigar.addWithCount(.mismatch, 1);
    try cigar.addWithCount(.match, 9);
    try cigar.addWithCount(.insertion, 1);

    try hits.list.append(try SearchHit(DNA).init(allocator, target1, cigar));

    var buffer: [1000]u8 = undefined;
    var fbs = std.io.fixedBufferStream(&buffer);
    try AlnoutWriter(DNA).write(fbs.writer(), query, hits);

    var it = std.mem.split(u8, fbs.getWritten(), "\n");
    try std.testing.expectEqualStrings("Query >QryId", it.next().?);
    try std.testing.expectEqualStrings(" %Id   TLen  Target", it.next().?);
    try std.testing.expectEqualStrings("100%     25  DbId3", it.next().?);
    try std.testing.expectEqualStrings(" 81%     33  DbId1", it.next().?);
    try std.testing.expectEqualStrings("", it.next().?);
    try std.testing.expectEqualStrings(" Query 27nt >QryId", it.next().?);
    try std.testing.expectEqualStrings("Target 25nt >DbId3", it.next().?);
    try std.testing.expectEqualStrings("", it.next().?);
    try std.testing.expectEqualStrings("Qry  1 + GGTGAGACGTTACGCAATAAATTGA 25", it.next().?);
    try std.testing.expectEqualStrings("         |||||||||||||||||||||||||   ", it.next().?);
    try std.testing.expectEqualStrings("Tgt  1 + GGTGAGACGTTACGCAATAAATTGA 25", it.next().?);

    // Query >QryId
    //  %Id   TLen  Target
    // 100%     25  DbId3
    //  81%     33  DbId1

    //  Query 27nt >QryId
    // Target 25nt >DbId3

    // Qry  1 + GGTGAGACGTTACGCAATAAATTGA 25
    //          |||||||||||||||||||||||||
    // Tgt  1 + GGTGAGACGTTACGCAATAAATTGA 25

    // 25 cols, 25 ids (100.0%), 0 gaps (0.0%)

    //  Query 27nt >QryId
    // Target 33nt >DbId1

    // Qry  1 + GGTGAGACGTTACGCAATAAATTGAGA 27
    //           ||||||||  | |||| |||||||||
    // Tgt  3 + CGTGAGACG--ATGCAAAAAATTGAGA 27

    // 27 cols, 22 ids (81.5%), 2 gaps (7.4%)

}