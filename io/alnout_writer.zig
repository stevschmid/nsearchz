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
    const MaxColumnsPerAlignmentLine = 60;

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

                var cigar = try hit.cigar.clone();
                defer cigar.deinit();

                const first = cigar.entries.items[0];
                switch (first.op) {
                    .deletion => {
                        target_idx = first.count;
                        _ = cigar.entries.orderedRemove(0);
                    },
                    .insertion => {
                        query_idx = first.count;
                        _ = cigar.entries.orderedRemove(0);
                    },
                    else => {},
                }

                const last = cigar.entries.items[cigar.entries.items.len - 1];
                switch (last.op) {
                    .deletion, .insertion => _ = cigar.entries.orderedRemove(cigar.entries.items.len - 1),
                    else => {},
                }

                const max_length = std.math.max(query.length(), hit.target.length());

                var top_cigar = try cigar.clone();
                defer top_cigar.deinit();
                var middle_cigar = try cigar.clone();
                defer middle_cigar.deinit();
                var bottom_cigar = try cigar.clone();
                defer bottom_cigar.deinit();

                while (!top_cigar.isEmpty()) {
                    // Qry  1 + GGTGAGACGTTACGCAATAAATTGA 25
                    try std.fmt.format(writer, "Qry ", .{});
                    try printLength(writer, max_length, query_idx + 1);

                    try std.fmt.format(writer, " + ", .{});

                    var count: usize = 0;
                    while (top_cigar.popOp()) |op| {
                        var ch = query.data[query_idx];
                        switch (op) {
                            .match, .mismatch, .insertion => query_idx += 1,
                            .deletion => ch = '-',
                        }

                        try std.fmt.format(writer, "{c}", .{ch});

                        count += 1;
                        if (count % MaxColumnsPerAlignmentLine == 0)
                            break;
                    }

                    try std.fmt.format(writer, " ", .{});
                    try printLength(writer, max_length, query_idx);
                    try std.fmt.format(writer, "\n", .{});

                    //          |||||||||||||||||||||||||
                    try std.fmt.format(writer, "    ", .{});
                    try printPadding(writer, max_length);
                    try std.fmt.format(writer, "   ", .{});

                    count = 0;
                    while (middle_cigar.popOp()) |op| {
                        const ch: u8 = switch (op) {
                            .match => '|',
                            else => ' ',
                        };
                        try std.fmt.format(writer, "{c}", .{ch});

                        count += 1;
                        if (count % MaxColumnsPerAlignmentLine == 0)
                            break;
                    }

                    try std.fmt.format(writer, " ", .{});
                    try printPadding(writer, max_length);
                    try std.fmt.format(writer, "\n", .{});

                    // Tgt  1 + GGTGAGACGTTACGCAATAAATTGA 25
                    try std.fmt.format(writer, "Tgt ", .{});
                    try printLength(writer, max_length, target_idx + 1);

                    try std.fmt.format(writer, " + ", .{});

                    count = 0;
                    while (bottom_cigar.popOp()) |op| {
                        var ch = hit.target.data[target_idx];
                        switch (op) {
                            .match, .mismatch, .deletion => target_idx += 1,
                            .insertion => ch = '-',
                        }

                        try std.fmt.format(writer, "{c}", .{ch});

                        count += 1;
                        if (count % MaxColumnsPerAlignmentLine == 0)
                            break;
                    }

                    try std.fmt.format(writer, " ", .{});
                    try printLength(writer, max_length, target_idx);
                    try std.fmt.format(writer, "\n", .{});

                    try std.fmt.format(writer, "\n", .{});
                }

                // Calc stats
                var num_cols: usize = 0;
                var num_matches: usize = 0;
                var num_gaps: usize = 0;

                while (cigar.popOp()) |op| {
                    switch (op) {
                        .match => num_matches += 1,
                        .insertion, .deletion => num_gaps += 1,
                        .mismatch => {},
                    }

                    num_cols += 1;
                }

                // 25 cols, 25 ids (100.0%), 0 gaps (0.0%)
                const ids_percent = 100.0 * @intToFloat(f32, num_matches) / @intToFloat(f32, num_cols);
                const gaps_percent = 100.0 * @intToFloat(f32, num_gaps) / @intToFloat(f32, num_cols);
                try std.fmt.format(writer, "{d} cols, {d} ids ({d:.1}%), {d} gaps ({d:.1}%)\n", .{ num_cols, num_matches, ids_percent, num_gaps, gaps_percent });

                try std.fmt.format(writer, "\n", .{});
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

const DNA = alphabet.DNA;

test "multiple hits per query" {
    const allocator = std.testing.allocator;

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
    try std.testing.expectEqualStrings("", it.next().?);
    try std.testing.expectEqualStrings("25 cols, 25 ids (100.0%), 0 gaps (0.0%)", it.next().?);
    try std.testing.expectEqualStrings("", it.next().?);
    try std.testing.expectEqualStrings(" Query 27nt >QryId", it.next().?);
    try std.testing.expectEqualStrings("Target 33nt >DbId1", it.next().?);
    try std.testing.expectEqualStrings("", it.next().?);
    try std.testing.expectEqualStrings("Qry  1 + GGTGAGACGTTACGCAATAAATTGAGA 27", it.next().?);
    try std.testing.expectEqualStrings("          ||||||||  | |||| |||||||||   ", it.next().?);
    try std.testing.expectEqualStrings("Tgt  3 + CGTGAGACG--ATGCAAAAAATTGAGA 27", it.next().?);
    try std.testing.expectEqualStrings("", it.next().?);
    try std.testing.expectEqualStrings("27 cols, 22 ids (81.5%), 2 gaps (7.4%)", it.next().?);
}

test "dna" {
    const allocator = std.testing.allocator;

    var cigar = Cigar.init(allocator);
    defer cigar.deinit();

    var hits = SearchHitList(DNA).init(allocator);
    defer hits.deinit();

    // 2D1X8=2I1=1X4=1X9=1I
    var query = try Sequence(DNA).init(allocator, "Procavia capensis", "CUUUGCCUGAACGCAAGACUCUUCAACCUCAGGACUUGCAGAAUUGGUAGAAUGCCGUCCUAAGGUUGUUGAGUUCUGUGUUUGGAGGC");
    defer query.deinit();

    var target = try Sequence(DNA).init(allocator, "Dasypus novemcinctus", "CGUCACCUGAACUCAUGACUCUUCAACUUCAGGACUUGCAGAAUUAAUGGAAUGCCGUCCUAAGGUUGUUGAGUUCUGCGUUUCUGGGC");
    defer target.deinit();

    try cigar.addFromString("1=1X1=2X7=1X2=1X11=1X17=2X1=1X29=1X4=3X3=");
    try hits.list.append(try SearchHit(DNA).init(allocator, target, cigar));

    var buffer: [1000]u8 = undefined;
    var fbs = std.io.fixedBufferStream(&buffer);
    try AlnoutWriter(DNA).write(fbs.writer(), query, hits);

    const expected =
        \\Query >Procavia capensis
        \\ %Id   TLen  Target
        \\ 85%     89  Dasypus novemcinctus
        \\
        \\ Query 89nt >Procavia capensis
        \\Target 89nt >Dasypus novemcinctus
        \\
        \\Qry  1 + CUUUGCCUGAACGCAAGACUCUUCAACCUCAGGACUUGCAGAAUUGGUAGAAUGCCGUCC 60
        \\         | |  ||||||| || ||||||||||| |||||||||||||||||  | |||||||||||   
        \\Tgt  1 + CGUCACCUGAACUCAUGACUCUUCAACUUCAGGACUUGCAGAAUUAAUGGAAUGCCGUCC 60
        \\
        \\Qry 61 + UAAGGUUGUUGAGUUCUGUGUUUGGAGGC 89
        \\         |||||||||||||||||| ||||   |||   
        \\Tgt 61 + UAAGGUUGUUGAGUUCUGCGUUUCUGGGC 89
        \\
        \\89 cols, 76 ids (85.4%), 0 gaps (0.0%)
        \\
        \\
    ;

    try std.testing.expectEqualStrings(expected, fbs.getWritten());
}

// auto entry = std::make_pair(
//   Sequence< DNA >( "RF00966;mir-676;ABRQ01840532.1/340-428   9813:Procavia "
//                    "capensis (cape rock hyrax)",
//                    "CUUUGCCUGAACGCAAGACUCUUCAACCUCAGGACUUGCAGAAUUGGUAGA"
//                    "AUGCCGUCCUAAGGUUGUUGAGUUCUGUGUUUGGAGGC" ),
//   HitList< DNA >( {
//     { { "RF00966;mir-676;AAGV020395671.1/1356-1444   9361:Dasypus "
//         "novemcinctus (nine-banded armadillo)",
//         "CGUCACCUGAACUCAUGACUCUUCAACUUCAGGACUUGCAGAAUUAAUGGAAUGCCGUCCUAAGGU"
//         "UGUUGAGUUCUGCGUUUCUGGGC" },
//       "1=1X1=2X7=1X2=1X11=1X17=2X1=1X29=1X4=3X3=" },
//   } ) );

// std::ostringstream oss;
// Alnout::Writer< DNA > writer( oss );

// SECTION( "Default" ) {
//   writer << entry;
//   REQUIRE( oss.str() == AlnoutOutputForDNA );
// }

// const char AlnoutOutputForDNA[] = R"(Query >RF00966;mir-676;ABRQ01840532.1/340-428   9813:Procavia capensis (cape rock hyrax)
//  %Id   TLen  Target
//  85%     89  RF00966;mir-676;AAGV020395671.1/1356-1444   9361:Dasypus novemcinctus (nine-banded armadillo)

//  Query 89nt >RF00966;mir-676;ABRQ01840532.1/340-428   9813:Procavia capensis (cape rock hyrax)
// Target 89nt >RF00966;mir-676;AAGV020395671.1/1356-1444   9361:Dasypus novemcinctus (nine-banded armadillo)

// Qry  1 + CUUUGCCUGAACGCAAGACUCUUCAACCUCAGGACUUGCAGAAUUGGUAGAAUGCCGUCC 60
//          | |  ||||||| || ||||||||||| |||||||||||||||||  | |||||||||||
// Tgt  1 + CGUCACCUGAACUCAUGACUCUUCAACUUCAGGACUUGCAGAAUUAAUGGAAUGCCGUCC 60

// Qry 61 + UAAGGUUGUUGAGUUCUGUGUUUGGAGGC 89
//          |||||||||||||||||| ||||   |||
// Tgt 61 + UAAGGUUGUUGAGUUCUGCGUUUCUGGGC 89

// 89 cols, 76 ids (85.4%), 0 gaps (0.0%)
