const std = @import("std");

const alphabet = @import("bio/alphabet.zig");
const Sequence = @import("sequence.zig").Sequence;
const Cigar = @import("cigar.zig").Cigar;
const CigarOp = @import("cigar.zig").CigarOp;

pub const BandedAlignOptions = struct {
    bandwidth: usize = 16,

    gap_interior_open_score: i32 = -20,
    gap_interior_extend_score: i32 = -2,

    gap_terminal_open_score: i32 = -2,
    gap_terminal_extend_score: i32 = -1,
};

pub const BandedAlignDirection = enum {
    forward,
    backward,
};

pub fn BandedAlign(comptime A: type) type {
    return struct {
        const Self = @This();
        const MinScore = -1_000_000;

        const Gap = struct {
            score: i32 = MinScore,
            is_terminal: bool = false,

            pub fn openOrExtend(self: *Gap, banded_align: *Self, score: i32, is_terminal: bool, length: usize) void {
                var new_gap_score: i32 = score;

                const o = banded_align.options;
                const gap_open_score = if (is_terminal) o.gap_terminal_open_score else o.gap_interior_open_score;
                const gap_extend_score = if (is_terminal) o.gap_terminal_extend_score else o.gap_interior_extend_score;

                if (length > 0) {
                    new_gap_score += gap_open_score + @intCast(i32, length) * gap_extend_score;
                }

                self.score += @intCast(i32, length) * gap_extend_score;

                if (new_gap_score > self.score) {
                    self.score = new_gap_score;
                    self.is_terminal = is_terminal;
                }
            }

            pub fn reset(self: *Gap) void {
                self.score = MinScore;
                self.is_terminal = false;
            }
        };

        options: BandedAlignOptions,
        allocator: std.mem.Allocator,
        scores: std.ArrayList(i32),
        vertical_gaps: std.ArrayList(Gap),
        ops: std.ArrayList(CigarOp),

        pub fn init(allocator: std.mem.Allocator, options: BandedAlignOptions) Self {
            return Self{
                .allocator = allocator,
                .options = options,
                .scores = std.ArrayList(i32).init(allocator),
                .vertical_gaps = std.ArrayList(Gap).init(allocator),
                .ops = std.ArrayList(CigarOp).init(allocator),
            };
        }

        pub fn deinit(self: *Self) void {
            self.scores.deinit();
            self.vertical_gaps.deinit();
            self.ops.deinit();
        }

        pub fn process(
            self: *Self,
            seq_one: Sequence(A),
            seq_two: Sequence(A),
            dir: BandedAlignDirection,
            start_one_: usize,
            start_two_: usize,
            end_one_: ?usize,
            end_two_: ?usize,
            cigar: ?*Cigar,
        ) !i32 {
            // const width = if (dir == AlignDirection.forward) (seq_one.data.len - start_one + 1) else (start_one + 1);
            // const height = if (dir == AlignDirection.forward) (seq_two.data.len - start_two + 1) else (start_two + 1);

            // try self.row.resize(@floatToInt(usize, @intToFloat(f32, width) * 1.5));
            // const row = self.row.items;

            // Calculate matrix width, depending on alignment
            // direction and length of sequences
            // A will be on the X axis (width of matrix)
            // B will be on the Y axis (height of matrix)
            var len_one = seq_one.data.len;
            var len_two = seq_two.data.len;

            const default_end_one = if (dir == BandedAlignDirection.forward) len_one else 0;
            var end_one = end_one_ orelse default_end_one;

            const default_end_two = if (dir == BandedAlignDirection.forward) len_two else 0;
            var end_two = end_two_ orelse default_end_two;

            var start_one = start_one_;
            var start_two = start_two_;

            // sanitize
            if (start_one > len_one) start_one = len_one;
            if (start_two > len_two) start_two = len_two;
            if (end_one > len_one) end_one = len_one;
            if (end_two > len_two) end_two = len_two;

            const width = if (end_one > start_one) end_one - start_one + 1 else start_one - end_one + 1;
            const height = if (end_two > start_two) end_two - start_two + 1 else start_two - end_two + 1;

            // Make sure we have enough cells
            try self.scores.resize(@floatToInt(usize, @intToFloat(f32, width) * 1.5));
            std.mem.set(i32, self.scores.items, MinScore);
            const scores = self.scores.items;

            try self.vertical_gaps.resize(@floatToInt(usize, @intToFloat(f32, width) * 1.5));
            std.mem.set(Gap, self.vertical_gaps.items, Gap{});
            const vertical_gaps = self.vertical_gaps.items;

            try self.ops.resize(@floatToInt(usize, @intToFloat(f32, width * height) * 1.5));
            const ops = self.ops.items;

            // Initialize first row
            const bw = self.options.bandwidth;

            const is_beginning_one: bool = (start_one == 0 or start_one == len_one);
            const is_beginning_two: bool = (start_two == 0 or start_two == len_two);

            const is_ending_one: bool = (end_one == 0 or end_one == len_one);
            const is_ending_two: bool = (end_two == 0 or end_two == len_two);

            scores[0] = 0;
            vertical_gaps[0].reset();
            vertical_gaps[0].openOrExtend(self, scores[0], is_beginning_two, 1);

            var horizontal_gap = Gap{};

            var x: usize = 1;
            while (x < width) : (x += 1) {
                if (x > bw and height > 1) // only break on BW bound if B is not empty
                    break;

                horizontal_gap.openOrExtend(self, scores[x - 1], is_beginning_one, 1);
                scores[x] = horizontal_gap.score;
                ops[x] = CigarOp.insertion;
                vertical_gaps[x].reset();
            }
            if (x < width) {
                scores[x] = MinScore;
                vertical_gaps[x].reset();
            }

            // Row by row
            var center: usize = 1;

            var y: usize = 1;
            while (y < height) : (y += 1) {
                var score: i32 = MinScore;

                // Calculate band bounds
                const left_bound = std.math.min(width - 1, if (center > bw) center - bw else 0);
                const right_bound = std.math.min(width - 1, center + bw);

                // Set diagonal score for first calculated cell in row
                var diag_score: i32 = MinScore;
                if (left_bound > 0) {
                    diag_score = scores[left_bound - 1];
                    scores[left_bound - 1] = MinScore;
                    vertical_gaps[left_bound - 1].reset();
                }

                // Calculate row within the band bounds
                horizontal_gap.reset();
                x = left_bound;
                while (x <= right_bound) : (x += 1) {
                    var pos_one: usize = 0;
                    var pos_two: usize = 0;

                    var match: bool = false;

                    if (x > 0) {
                        // diagScore: score at col-1, row-1
                        pos_one = if (dir == BandedAlignDirection.forward) start_one + x - 1 else start_one - x;
                        pos_two = if (dir == BandedAlignDirection.forward) start_two + y - 1 else start_two - y;

                        const letter_one = seq_one.data[pos_one];
                        const letter_two = seq_two.data[pos_two];

                        match = A.match(letter_one, letter_two);
                        score = diag_score + A.score(letter_one, letter_two);
                    }

                    // Select highest score
                    //  - coming from diag (current),
                    //  - coming from left (row)
                    //  - coming from top (col)
                    const vertical_gap = &vertical_gaps[x];
                    score = std.math.max(score, horizontal_gap.score);
                    score = std.math.max(score, vertical_gap.score);

                    // Save the prev score at (x - 1) which
                    // we will use to compute the diagonal score at (x)
                    diag_score = scores[x];

                    // Save new score
                    scores[x] = score;

                    // Record op
                    var op: CigarOp = undefined;
                    if (score == horizontal_gap.score) {
                        op = CigarOp.insertion;
                    } else if (score == vertical_gap.score) {
                        op = CigarOp.deletion;
                    } else if (match) {
                        op = CigarOp.match;
                    } else {
                        op = CigarOp.mismatch;
                    }
                    ops[y * width + x] = op;

                    // Calculate potential gaps
                    const is_gap_terminal_one = (x == 0 or x == width - 1) and is_ending_one;
                    const is_gap_terminal_two = (y == 0 or y == height - 1) and is_ending_two;

                    horizontal_gap.openOrExtend(self, score, is_gap_terminal_two, 1);
                    vertical_gap.openOrExtend(self, score, is_gap_terminal_one, 1);
                }

                if (right_bound + 1 < width) {
                    scores[right_bound + 1] = MinScore;
                    vertical_gaps[right_bound + 1].reset();
                }

                if (right_bound == left_bound)
                    break;

                // Move one cell over for the next row
                center += 1;
            }

            // backtrack
            if (cigar != null) {
                var bx = x - 1;
                var by = y - 1;

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

                cigar.?.reverse();
            }

            var score = scores[x - 1];
            if (x == width) {
                // We reached the end of A, emulate going down on B (vertical gaps)
                const remaining_two = height - y;
                const vertical_gap = &vertical_gaps[x - 1];
                vertical_gap.openOrExtend(self, score, vertical_gap.is_terminal, remaining_two);
                score = vertical_gap.score;
            } else if (y == height) {
                // We reached the end of B, emulate going down on A (horizontal gaps)
                const remaining_one = width - x;
                horizontal_gap.openOrExtend(self, score, horizontal_gap.is_terminal, remaining_one);
                score = horizontal_gap.score;
            }

            _ = cigar;

            return score;
        }
    };
}

test "basic" {
    const allocator = std.testing.allocator;

    var seq_one = try Sequence(alphabet.DNA).init(allocator, "one", "TATAATGTTTACATTGG");
    defer seq_one.deinit();

    var seq_two = try Sequence(alphabet.DNA).init(allocator, "two", "TATAATGACACTGG");
    defer seq_two.deinit();

    const options = BandedAlignOptions{};
    var banded_align = BandedAlign(alphabet.DNA).init(allocator, options);
    defer banded_align.deinit();

    var cigar = Cigar.init(allocator);
    defer cigar.deinit();

    var score = try banded_align.process(seq_one, seq_two, BandedAlignDirection.forward, 0, 0, null, null, &cigar);

    // TATAATGTTTACATTGG
    // |||||||   |||.|||
    // TATAATG---ACACTGG
    const MatchScore = 2;
    const MismatchScore = -4;

    try std.testing.expectEqual(13 * MatchScore + 1 * options.gap_interior_open_score + 3 * options.gap_interior_extend_score + 1 * MismatchScore, score);

    var cigar_str = try cigar.toStringAlloc(allocator);
    defer allocator.free(cigar_str);
    try std.testing.expectEqualStrings("7=3I3=1X3=", cigar_str);
}

test "gap penalties" {
    const allocator = std.testing.allocator;

    var seq_one = try Sequence(alphabet.DNA).init(allocator, "one", "GGATCCTA");
    defer seq_one.deinit();

    var seq_two = try Sequence(alphabet.DNA).init(allocator, "two", "ATCGTA");
    defer seq_two.deinit();

    {
        var banded_align = BandedAlign(alphabet.DNA).init(allocator, BandedAlignOptions{});
        defer banded_align.deinit();

        var cigar = Cigar.init(allocator);
        defer cigar.deinit();

        _ = try banded_align.process(seq_one, seq_two, BandedAlignDirection.forward, 0, 0, null, null, &cigar);
        var cigar_str = try cigar.toStringAlloc(allocator);
        defer allocator.free(cigar_str);
        try std.testing.expectEqualStrings("2I3=1X2=", cigar_str);
    }

    {
        var banded_align = BandedAlign(alphabet.DNA).init(allocator, BandedAlignOptions{ .gap_terminal_open_score = -40 });
        defer banded_align.deinit();

        var cigar = Cigar.init(allocator);
        defer cigar.deinit();

        _ = try banded_align.process(seq_one, seq_two, BandedAlignDirection.forward, 0, 0, null, null, &cigar);
        var cigar_str = try cigar.toStringAlloc(allocator);
        defer allocator.free(cigar_str);
        try std.testing.expectEqualStrings("1X2I2=1X2=", cigar_str);
    }
}
