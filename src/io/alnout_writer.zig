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
        const MaxColumnsPerAlignmentLine = 60;

        pub fn write(unbuffered_writer: anytype, query: Sequence(A), hits: SearchHitList(A)) !void {
            // buffer, IO directly is slow
            var buffered_writer = std.io.bufferedWriter(unbuffered_writer);
            var writer = buffered_writer.writer();

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

                var top_iter = cigar.iterator();
                var middle_iter = cigar.iterator();
                var bottom_iter = cigar.iterator();

                while (!top_iter.isEmpty()) {
                    const query_idx_line_start = query_idx;
                    const target_idx_line_start = target_idx;

                    // Qry  1 + GGTGAGACGTTACGCAATAAATTGA 25
                    try std.fmt.format(writer, "Qry ", .{});
                    try printLength(writer, max_length, queryPos(query_idx, query, hit), .left);

                    if (A.SupportsStrands) {
                        try std.fmt.format(writer, " {c} ", .{queryStrand(hit)});
                    } else {
                        try std.fmt.format(writer, " ", .{});
                    }

                    var count: usize = 0;
                    while (top_iter.next()) |op| {
                        var ch = queryLetter(query_idx, query, hit);
                        switch (op) {
                            .match, .mismatch, .insertion => query_idx += 1,
                            .deletion => ch = '-',
                        }

                        try writer.writeByte(ch);

                        count += 1;
                        if (count % MaxColumnsPerAlignmentLine == 0)
                            break;
                    }

                    try std.fmt.format(writer, " ", .{});
                    try printLength(writer, max_length, queryPos(query_idx - 1, query, hit), .right);
                    try std.fmt.format(writer, "\n", .{});

                    //          |||||||||||||||||||||||||
                    try std.fmt.format(writer, "    ", .{});
                    try printPadding(writer, max_length);

                    if (A.SupportsStrands) {
                        try std.fmt.format(writer, "   ", .{});
                    } else {
                        try std.fmt.format(writer, " ", .{});
                    }

                    var query_match_idx = query_idx_line_start;
                    var target_match_idx = target_idx_line_start;

                    count = 0;
                    while (middle_iter.next()) |op| {
                        const ch: u8 = sym: {
                            switch (op) {
                                .match, .mismatch => {
                                    const query_letter = queryLetter(query_match_idx, query, hit);
                                    const target_letter = targetLetter(target_match_idx, hit);
                                    query_match_idx += 1;
                                    target_match_idx += 1;
                                    break :sym A.matchSymbol(query_letter, target_letter);
                                },
                                .deletion => {
                                    target_match_idx += 1;
                                    break :sym ' ';
                                },
                                .insertion => {
                                    query_match_idx += 1;
                                    break :sym ' ';
                                },
                            }
                        };

                        try writer.writeByte(ch);

                        count += 1;
                        if (count % MaxColumnsPerAlignmentLine == 0)
                            break;
                    }

                    try std.fmt.format(writer, " ", .{});
                    try printPadding(writer, max_length);
                    try std.fmt.format(writer, "\n", .{});

                    // Tgt  1 + GGTGAGACGTTACGCAATAAATTGA 25
                    try std.fmt.format(writer, "Tgt ", .{});
                    try printLength(writer, max_length, targetPos(target_idx, hit), .left);

                    if (A.SupportsStrands) {
                        try std.fmt.format(writer, " {c} ", .{targetStrand(hit)});
                    } else {
                        try std.fmt.format(writer, " ", .{});
                    }

                    count = 0;
                    while (bottom_iter.next()) |op| {
                        var ch = targetLetter(target_idx, hit);
                        switch (op) {
                            .match, .mismatch, .deletion => target_idx += 1,
                            .insertion => ch = '-',
                        }

                        try writer.writeByte(ch);

                        count += 1;
                        if (count % MaxColumnsPerAlignmentLine == 0)
                            break;
                    }

                    try std.fmt.format(writer, " ", .{});
                    try printLength(writer, max_length, targetPos(target_idx - 1, hit), .right);
                    try std.fmt.format(writer, "\n", .{});

                    try std.fmt.format(writer, "\n", .{});
                }

                // Calc stats
                var num_cols: usize = 0;
                var num_matches: usize = 0;
                var num_gaps: usize = 0;

                var iter = cigar.iterator();
                while (iter.next()) |op| {
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

                // flush content
                try buffered_writer.flush();
            }
        }

        fn queryPos(index: usize, query: Sequence(A), hit: SearchHit(A)) usize {
            if (hit.reverse_match) {
                return query.length() - index;
            } else {
                return index + 1;
            }
        }

        fn queryLetter(index: usize, query: Sequence(A), hit: SearchHit(A)) u8 {
            if (hit.reverse_match) {
                return A.complement(query.data[query.length() - index - 1]);
            } else {
                return query.data[index];
            }
        }

        fn queryStrand(hit: SearchHit(A)) u8 {
            return if (hit.reverse_match) '-' else '+';
        }

        fn targetPos(index: usize, hit: SearchHit(A)) usize {
            _ = hit;
            return index + 1;
        }

        fn targetLetter(index: usize, hit: SearchHit(A)) u8 {
            _ = hit;
            return hit.target.data[index];
        }

        fn targetStrand(hit: SearchHit(A)) u8 {
            _ = hit;
            return '+';
        }

        fn printPadding(writer: anytype, length: usize) !void {
            const max = std.fmt.count("{d}", .{length});
            var cur: usize = 0;

            while (cur < max) : (cur += 1) {
                try std.fmt.format(writer, " ", .{});
            }
        }

        const Pad = enum { left, right };
        fn printLength(writer: anytype, max_length: usize, length: usize, pad: Pad) !void {
            const max = std.fmt.count("{d}", .{max_length});
            var cur = std.fmt.count("{d}", .{length});

            if (pad == .left) try writer.writeByteNTimes(' ', max - cur);
            try std.fmt.format(writer, "{d}", .{length});
            if (pad == .right) try writer.writeByteNTimes(' ', max - cur);
        }
    };
}

const DNA = alphabet.DNA;
const Protein = alphabet.Protein;

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

    try hits.list.append(try SearchHit(DNA).init(allocator, target3, cigar, false));

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

    try hits.list.append(try SearchHit(DNA).init(allocator, target1, cigar, false));

    var buffer: [4096]u8 = undefined;
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
    try hits.list.append(try SearchHit(DNA).init(allocator, target, cigar, false));

    var buffer: [4096]u8 = undefined;
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

test "dna minus" {
    const allocator = std.testing.allocator;

    var cigar = Cigar.init(allocator);
    defer cigar.deinit();

    var hits = SearchHitList(DNA).init(allocator);
    defer hits.deinit();

    // 2D1X8=2I1=1X4=1X9=1I
    var query = try Sequence(DNA).init(allocator, "RevComp of Procavia capensis", "GCCTCCAAACACAGAACTCAACAACCTTAGGACGGCATTCTACCAATTCTGCAAGTCCTGAGGTTGAAGAGTCTTGCGTTCAGGCAAAG");
    defer query.deinit();

    var target = try Sequence(DNA).init(allocator, "Dasypus novemcinctus", "CGUCACCUGAACUCAUGACUCUUCAACUUCAGGACUUGCAGAAUUAAUGGAAUGCCGUCCUAAGGUUGUUGAGUUCUGCGUUUCUGGGC");
    defer target.deinit();

    try cigar.addFromString("1=1X1=2X7=1X2=1X11=1X17=2X1=1X29=1X4=3X3=");
    try hits.list.append(try SearchHit(DNA).init(allocator, target, cigar, true));

    var buffer: [4096]u8 = undefined;
    var fbs = std.io.fixedBufferStream(&buffer);
    try AlnoutWriter(DNA).write(fbs.writer(), query, hits);

    const expected =
        \\Query >RevComp of Procavia capensis
        \\ %Id   TLen  Target
        \\ 85%     89  Dasypus novemcinctus
        \\
        \\ Query 89nt >RevComp of Procavia capensis
        \\Target 89nt >Dasypus novemcinctus
        \\
        \\Qry 89 - CTTTGCCTGAACGCAAGACTCTTCAACCTCAGGACTTGCAGAATTGGTAGAATGCCGTCC 30
        \\         | +  ||+|||| || |||+|++|||| +||||||++||||||++  + |||+||||+||   
        \\Tgt  1 + CGUCACCUGAACUCAUGACUCUUCAACUUCAGGACUUGCAGAAUUAAUGGAAUGCCGUCC 60
        \\
        \\Qry 29 - TAAGGTTGTTGAGTTCTGTGTTTGGAGGC 1 
        \\         +||||++|++|||++|+| |+++   |||   
        \\Tgt 61 + UAAGGUUGUUGAGUUCUGCGUUUCUGGGC 89
        \\
        \\89 cols, 76 ids (85.4%), 0 gaps (0.0%)
        \\
        \\
    ;

    try std.testing.expectEqualStrings(expected, fbs.getWritten());
}

test "protein" {
    const allocator = std.testing.allocator;

    // auto entry1 = std::make_pair( Sequence< Protein >( "query1", "LAFQGVRN" ),
    //                               HitList< Protein >( {
    //                                 { { "target50", "MAFQGVRS" }, "1X6=1X" },
    //                                 { { "target114", "LAGQGSAN" }, "4=3X1=" },
    //                               } ) );
    // auto entry2 =
    //   std::make_pair( Sequence< Protein >( "query2", "GGGGGYFDEATGVCPF" ),
    //                   HitList< Protein >( {
    //                     { { "target1337", "YFDEATGICPFQQQ" }, "5I7=1X3=3D" },
    //                   } ) );

    // std::ostringstream oss;
    // Alnout::Writer< Protein > writer( oss );

    // writer << entry1;
    // writer << entry2;

    // REQUIRE( oss.str() == AlnoutOutputForProtein );

    var buffer: [4096]u8 = undefined;
    var fbs = std.io.fixedBufferStream(&buffer);

    {
        var hits = SearchHitList(Protein).init(allocator);
        defer hits.deinit();

        var query = try Sequence(Protein).init(allocator, "query1", "LAFQGVRN");
        defer query.deinit();

        var target50 = try Sequence(Protein).init(allocator, "target50", "MAFQGVRS");
        defer target50.deinit();
        var cigar50 = Cigar.init(allocator);
        defer cigar50.deinit();
        try cigar50.addFromString("1X6=1X");
        try hits.list.append(try SearchHit(Protein).init(allocator, target50, cigar50, false));

        var target114 = try Sequence(Protein).init(allocator, "target114", "LAGQGSAN");
        defer target114.deinit();
        var cigar114 = Cigar.init(allocator);
        defer cigar114.deinit();
        try cigar114.addFromString("4=3X1=");
        try hits.list.append(try SearchHit(Protein).init(allocator, target114, cigar114, false));

        try AlnoutWriter(Protein).write(fbs.writer(), query, hits);
    }

    {
        var hits = SearchHitList(Protein).init(allocator);
        defer hits.deinit();

        var query = try Sequence(Protein).init(allocator, "query2", "GGGGGYFDEATGVCPF");
        defer query.deinit();

        var target1337 = try Sequence(Protein).init(allocator, "target1337", "YFDEATGICPFQQQ");
        defer target1337.deinit();
        var cigar1337 = Cigar.init(allocator);
        defer cigar1337.deinit();
        try cigar1337.addFromString("5I7=1X3=3D");
        try hits.list.append(try SearchHit(Protein).init(allocator, target1337, cigar1337, false));

        try AlnoutWriter(Protein).write(fbs.writer(), query, hits);
    }

    const expected =
        \\Query >query1
        \\ %Id   TLen  Target
        \\ 75%      8  target50
        \\ 63%      8  target114
        \\
        \\ Query 8aa >query1
        \\Target 8aa >target50
        \\
        \\Qry 1 LAFQGVRN 8
        \\      :||||||.  
        \\Tgt 1 MAFQGVRS 8
        \\
        \\8 cols, 6 ids (75.0%), 0 gaps (0.0%)
        \\
        \\ Query 8aa >query1
        \\Target 8aa >target114
        \\
        \\Qry 1 LAFQGVRN 8
        \\      || ||  |  
        \\Tgt 1 LAGQGSAN 8
        \\
        \\8 cols, 5 ids (62.5%), 0 gaps (0.0%)
        \\
        \\Query >query2
        \\ %Id   TLen  Target
        \\ 91%     14  target1337
        \\
        \\ Query 16aa >query2
        \\Target 14aa >target1337
        \\
        \\Qry  6 YFDEATGVCPF 16
        \\       |||||||:|||   
        \\Tgt  1 YFDEATGICPF 11
        \\
        \\11 cols, 10 ids (90.9%), 0 gaps (0.0%)
        \\
        \\
    ;

    try std.testing.expectEqualStrings(expected, fbs.getWritten());
}
