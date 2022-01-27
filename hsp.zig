const std = @import("std");

pub const HSP = struct {
    const Self = @This();

    a_start: usize,
    a_end: usize,

    b_start: usize,
    b_end: usize,

    // cigar: Cigar,

    score: i32 = 0,

    pub fn length(self: Self) usize {
        return std.math.max(self.a_end - self.a_start, self.b_end - self.b_start) + 1;
    }

    pub fn is_overlapping(self: Self, other: HSP) bool {
        if (self.a_start <= other.a_end and other.a_start <= self.a_end) // overlap in A direction
            return true;

        if (self.b_start <= other.b_end and other.b_start <= self.b_end) // overlap in B direction
            return true;

        return false;
    }

    pub fn distance_to(self: Self, other: HSP) usize {
        var dx = if (self.a_start > other.a_end) self.a_start - other.a_end else other.a_start - self.a_end;
        if (dx > 0)
            dx -= 1;

        var dy = if (self.b_start > other.b_end) self.b_start - other.b_end else other.b_start - self.b_end;
        if (dy > 0)
            dy -= 1;

        return std.math.sqrt(dx*dx + dy*dy);
    }
};

test "length" {
    const hsp = HSP { .a_start = 5, .a_end = 6, .b_start = 11, .b_end = 15 };
    try std.testing.expectEqual(@as(usize, 5), hsp.length());

    const hsp2 = HSP { .a_start = 1, .a_end = 6, .b_start = 11, .b_end = 15 };
    try std.testing.expectEqual(@as(usize, 6), hsp2.length());
}

test "overlapping" {
    const hsp = HSP { .a_start = 1, .a_end = 2, .b_start = 25, .b_end = 27 };

    try std.testing.expectEqual(true, hsp.is_overlapping(HSP { .a_start = 1, .a_end = 2, .b_start = 55, .b_end = 57 } ));
    try std.testing.expectEqual(true, hsp.is_overlapping(HSP { .a_start = 2, .a_end = 3, .b_start = 55, .b_end = 57 } ));
    try std.testing.expectEqual(false, hsp.is_overlapping(HSP { .a_start = 3, .a_end = 4, .b_start = 55, .b_end = 57 } ));
    try std.testing.expectEqual(true, hsp.is_overlapping(HSP { .a_start = 3, .a_end = 4, .b_start = 20, .b_end = 28 } ));
    try std.testing.expectEqual(true, hsp.is_overlapping(HSP { .a_start = 3, .a_end = 4, .b_start = 20, .b_end = 25 } ));
    try std.testing.expectEqual(false, hsp.is_overlapping(HSP { .a_start = 3, .a_end = 4, .b_start = 20, .b_end = 24 } ));
    try std.testing.expectEqual(true, hsp.is_overlapping(HSP { .a_start = 0, .a_end = 100, .b_start = 0, .b_end = 100 }));
}

test "distance" {
    {
        const hsp1 = HSP { .a_start = 1, .a_end = 1, .b_start = 25, .b_end = 25 };
        const hsp2 = HSP { .a_start = 2, .a_end = 2, .b_start = 26, .b_end = 26 };
        try std.testing.expectEqual(@as(usize, 0), hsp1.distance_to(hsp2));
    }

    {
        const hsp1 = HSP { .a_start = 1, .a_end = 2, .b_start = 25, .b_end = 26 };
        const hsp2 = HSP { .a_start = 5, .a_end = 10, .b_start = 30, .b_end = 35 };
        // x: 2 -> 5 = 2 empty cells
        // y: 26 -> 30 = 3 empty cells
        // => sqrt(2*2 + 3*3)
        try std.testing.expectEqual(@as(usize, 3), hsp1.distance_to(hsp2));
        try std.testing.expectEqual(@as(usize, 3), hsp2.distance_to(hsp1));
    }
}
