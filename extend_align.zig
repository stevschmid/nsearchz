const std = @import("std");

const alphabet = @import("bio/alphabet.zig");
const Sequence = @import("sequence.zig").Sequence;
const Cigar = @import("cigar.zig").Cigar;
const CigarOp = @import("cigar.zig").CigarOp;

pub const ExtendAlignResult = struct {
    score: i32,
    pos_one: usize,
    pos_two: usize,
};

pub const ExtendAlignDirection = enum {
    forward,
    backward,
};

pub const ExtendAlignOptions = struct {
    x_drop: i32 = 32,
    gap_open_score: i32 = -20,
    gap_extend_score: i32 = -2,
};

pub fn ExtendAlign(comptime A: type) type {
    return struct {
        const Self = @This();
        const MinScore = -1_000_000;

        const Cell = struct {
            score: i32,
            score_gap: i32,
        };

        options: ExtendAlignOptions,
        allocator: std.mem.Allocator,
        row: std.ArrayList(Cell),
        ops: std.ArrayList(CigarOp),

        pub fn init(allocator: std.mem.Allocator, options: ExtendAlignOptions) Self {
            return Self{
                .allocator = allocator,
                .options = options,
                .row = std.ArrayList(Cell).init(allocator),
                .ops = std.ArrayList(CigarOp).init(allocator),
            };
        }

        pub fn deinit(self: *Self) void {
            self.row.deinit();
            self.ops.deinit();
        }

        pub fn extend(self: *Self, seq_one: Sequence(A), seq_two: Sequence(A), dir: ExtendAlignDirection, start_one: usize, start_two: usize, cigar: ?*Cigar) !ExtendAlignResult {
            const width = if (dir == ExtendAlignDirection.forward) (seq_one.data.len - start_one + 1) else (start_one + 1);
            const height = if (dir == ExtendAlignDirection.forward) (seq_two.data.len - start_two + 1) else (start_two + 1);

            try self.row.resize(@floatToInt(usize, @intToFloat(f32, width) * 1.5));
            const row = self.row.items;

            try self.ops.resize(@floatToInt(usize, @intToFloat(f32, width * height) * 1.5));
            const ops = self.ops.items;

            var x_drop = self.options.x_drop;
            var gap_open_score = self.options.gap_open_score;
            var gap_extend_score = self.options.gap_extend_score;

            var best_one: usize = start_one;
            var best_two: usize = start_two;

            var best_x: usize = 0;
            var best_y: usize = 0;
            var best_score: i32 = 0;

            // init row
            row[0].score = 0;
            row[0].score_gap = gap_open_score + gap_extend_score;

            var score: i32 = 0;
            var x: usize = 1;
            while (x < width) : (x += 1) {
                score = gap_open_score + @intCast(i32, x) * gap_extend_score;

                if (score < -x_drop)
                    break;

                ops[x] = CigarOp.insertion;
                row[x].score = score;
                row[x].score_gap = MinScore;
            }

            var row_size: usize = x;
            var first_x: usize = 0;

            // row by row
            var y: usize = 1;
            while (y < height) : (y += 1) {
                var row_gap: i32 = MinScore;
                var diag_score: i32 = MinScore;
                score = MinScore;

                var last_x: usize = first_x;

                x = first_x;
                while (x < row_size) : (x += 1) {
                    var col_gap = row[x].score_gap;
                    _ = col_gap;

                    var pos_one: usize = 0;
                    var pos_two: usize = 0;

                    var match: bool = false;

                    if (x > 0) {
                        // diagScore: score at col-1, row-1
                        pos_one = if (dir == ExtendAlignDirection.forward) start_one + x - 1 else start_one - x;
                        pos_two = if (dir == ExtendAlignDirection.forward) start_two + y - 1 else start_two - y;

                        const letter_one = seq_one.data[pos_one];
                        const letter_two = seq_two.data[pos_two];

                        match = A.match(letter_one, letter_two);
                        score = diag_score + A.score(letter_one, letter_two);
                    }

                    // select highest score
                    //  - coming from diag (current),
                    //  - coming from left (row)
                    //  - coming from top (col)
                    score = std.math.max(score, row_gap);
                    score = std.math.max(score, col_gap);

                    // row[x] right now points to the previous row, so use this
                    // in the next iteration for the diagonal computation of (x, y )
                    diag_score = row[x].score;

                    if (best_score - score > x_drop) {
                        // X-Drop test failed
                        row[x].score = MinScore;

                        if (x == first_x) {
                            // Tighten left bound
                            first_x += 1;
                        }
                    } else {
                        last_x = x;

                        // Check if we achieved new highscore
                        if (score > best_score) {
                            best_score = score;

                            best_x = x;
                            best_y = y;

                            best_one = pos_one;
                            best_two = pos_two;
                        }

                        // Record new score
                        var op: CigarOp = undefined;
                        if (score == row_gap) {
                            op = CigarOp.insertion;
                        } else if (score == col_gap) {
                            op = CigarOp.deletion;
                        } else if (match) {
                            op = CigarOp.match;
                        } else {
                            op = CigarOp.mismatch;
                        }
                        ops[y * width + x] = op;

                        // update scores
                        row[x].score = score;
                        row[x].score_gap = std.math.max(score + gap_open_score + gap_extend_score, col_gap + gap_extend_score);
                        row_gap = std.math.max(score + gap_open_score + gap_extend_score, row_gap + gap_extend_score);
                    }
                } // while x

                if (first_x == row_size) {
                    // All cells failed the X-Drop test
                    // We are done 8)
                    break;
                }

                if (last_x < row_size - 1) {
                    // Tighten right bound
                    row_size = last_x + 1;
                } else {
                    // Extend row, since last checked column didn't fail X-Drop test
                    while (row_gap >= (best_score - x_drop) and row_size < width) {
                        row[row_size].score = row_gap;
                        row[row_size].score_gap = row_gap + gap_open_score + gap_extend_score;
                        ops[y * width + row_size] = CigarOp.insertion;
                        row_gap += gap_extend_score;
                        row_size += 1;
                    }
                }

                // Properly reset right bound
                if (row_size < width) {
                    row[row_size].score = MinScore;
                    row[row_size].score_gap = MinScore;
                    row_size += 1;
                }
            } // while y

            // backtrack
            if (cigar != null) {
                var bx = best_x;
                var by = best_y;

                while (bx > 0 or by > 0) {
                    const op = ops[by * width + bx];
                    try cigar.?.add(op);

                    // where did we come from?
                    switch (op) {
                        CigarOp.insertion => {
                            bx -= 1;
                        },
                        CigarOp.deletion => {
                            by -= 1;
                        },
                        CigarOp.match, CigarOp.mismatch => {
                            bx -= 1;
                            by -= 1;
                        },
                    }
                }

                if (dir == ExtendAlignDirection.forward) {
                    cigar.?.reverse();
                }
            }

            return ExtendAlignResult{
                .score = best_score,
                .pos_one = best_one,
                .pos_two = best_two,
            };
        }
    };
}

test "forward" {
    const allocator = std.testing.allocator;

    var seq_one = try Sequence(alphabet.DNA).init(allocator, "one", "GATTGCGGGG");
    defer seq_one.deinit();

    var seq_two = try Sequence(alphabet.DNA).init(allocator, "two", "GAGCGGT");
    defer seq_two.deinit();

    var extend_align = ExtendAlign(alphabet.DNA).init(allocator, ExtendAlignOptions{});
    defer extend_align.deinit();

    var cigar = Cigar.init(allocator);
    defer cigar.deinit();

    var result = try extend_align.extend(seq_one, seq_two, ExtendAlignDirection.forward, 0, 0, &cigar);
    try std.testing.expectEqual(@as(i32, 4), result.score); // 2 matches = +4
    // try std.testing.expectEqualSlices(CigarOp, &[_]CigarOp{ CigarOp.match, CigarOp.match }, cigar.ops.items);

}

test "gaps" {
    const allocator = std.testing.allocator;

    var seq_one = try Sequence(alphabet.DNA).init(allocator, "one", "GATTGCGGGG");
    defer seq_one.deinit();

    var seq_two = try Sequence(alphabet.DNA).init(allocator, "two", "GAGCGGT");
    defer seq_two.deinit();

    var extend_align = ExtendAlign(alphabet.DNA).init(allocator, ExtendAlignOptions{ .gap_open_score = -3 });
    defer extend_align.deinit();

    var cigar = Cigar.init(allocator);
    defer cigar.deinit();

    var result = try extend_align.extend(seq_one, seq_two, ExtendAlignDirection.forward, 0, 0, &cigar);
    try std.testing.expectEqual(@as(i32, 5), result.score); // 6 matches, 1 gap,  gap len = 1 - 3 - *

    var cigar_str = try cigar.toStringAlloc(allocator);
    defer allocator.free(cigar_str);
    try std.testing.expectEqualStrings("2=2I4=", cigar_str);
}

test "forward extend" {
    const allocator = std.testing.allocator;

    var seq_one = try Sequence(alphabet.DNA).init(allocator, "one", "ATCGG");
    defer seq_one.deinit();

    var seq_two = try Sequence(alphabet.DNA).init(allocator, "two", "ATCGT");
    defer seq_two.deinit();

    var extend_align = ExtendAlign(alphabet.DNA).init(allocator, ExtendAlignOptions{});
    defer extend_align.deinit();

    {
        var cigar = Cigar.init(allocator);
        defer cigar.deinit();

        var result = try extend_align.extend(seq_one, seq_two, ExtendAlignDirection.forward, 0, 0, &cigar);
        try std.testing.expectEqual(result.pos_one, 3);
        try std.testing.expectEqual(result.pos_two, 3);

        var cigar_str = try cigar.toStringAlloc(allocator);
        defer allocator.free(cigar_str);

        try std.testing.expectEqualStrings("4=", cigar_str);
    }

    {
        var cigar = Cigar.init(allocator);
        defer cigar.deinit();

        var result = try extend_align.extend(seq_one, seq_two, ExtendAlignDirection.forward, 3, 3, &cigar);
        try std.testing.expectEqual(result.pos_one, 3);
        try std.testing.expectEqual(result.pos_two, 3);

        var cigar_str = try cigar.toStringAlloc(allocator);
        defer allocator.free(cigar_str);

        try std.testing.expectEqualStrings("1=", cigar_str);
    }

    {
        var cigar = Cigar.init(allocator);
        defer cigar.deinit();

        var result = try extend_align.extend(seq_one, seq_two, ExtendAlignDirection.forward, 4, 4, &cigar);
        try std.testing.expectEqual(result.pos_one, 4);
        try std.testing.expectEqual(result.pos_two, 4);

        var cigar_str = try cigar.toStringAlloc(allocator);
        defer allocator.free(cigar_str);

        try std.testing.expectEqualStrings("", cigar_str);
    }
}

test "backward" {
    const allocator = std.testing.allocator;

    var seq_one = try Sequence(alphabet.DNA).init(allocator, "one", "ATCGGTTG");
    defer seq_one.deinit();

    var seq_two = try Sequence(alphabet.DNA).init(allocator, "two", "TCGGTAT");
    defer seq_two.deinit();

    var extend_align = ExtendAlign(alphabet.DNA).init(allocator, ExtendAlignOptions{});
    defer extend_align.deinit();

    var result = try extend_align.extend(seq_one, seq_two, ExtendAlignDirection.backward, 3, 2, null);
    try std.testing.expectEqual(@as(usize, 1), result.pos_one);
    try std.testing.expectEqual(@as(usize, 0), result.pos_two);
}
