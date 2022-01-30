const std = @import("std");

const Cigar = @import("cigar.zig").Cigar;

pub const HSP = struct {
    const Self = @This();

    start_one: usize,
    end_one: usize,

    start_two: usize,
    end_two: usize,

    pub fn length(self: Self) usize {
        return std.math.max(self.end_one - self.start_one, self.end_two - self.start_two) + 1;
    }

    pub fn is_overlapping(self: Self, other: HSP) bool {
        if (self.start_one <= other.end_one and other.start_one <= self.end_one) // overlap in A direction
            return true;

        if (self.start_two <= other.end_two and other.start_two <= self.end_two) // overlap in B direction
            return true;

        return false;
    }

    pub fn distance_to(self: Self, other: HSP) usize {
        var dx = if (self.start_one > other.end_one) self.start_one - other.end_one else other.start_one - self.end_one;
        if (dx > 0)
            dx -= 1;

        var dy = if (self.start_two > other.end_two) self.start_two - other.end_two else other.start_two - self.end_two;
        if (dy > 0)
            dy -= 1;

        return std.math.sqrt(dx * dx + dy * dy);
    }
};

test "length" {
    const hsp = HSP{ .start_one = 5, .end_one = 6, .start_two = 11, .end_two = 15 };
    try std.testing.expectEqual(@as(usize, 5), hsp.length());

    const hsp2 = HSP{ .start_one = 1, .end_one = 6, .start_two = 11, .end_two = 15 };
    try std.testing.expectEqual(@as(usize, 6), hsp2.length());
}

test "overlapping" {
    const hsp = HSP{ .start_one = 1, .end_one = 2, .start_two = 25, .end_two = 27 };

    try std.testing.expectEqual(true, hsp.is_overlapping(.{ .start_one = 1, .end_one = 2, .start_two = 55, .end_two = 57 }));
    try std.testing.expectEqual(true, hsp.is_overlapping(.{ .start_one = 2, .end_one = 3, .start_two = 55, .end_two = 57 }));
    try std.testing.expectEqual(false, hsp.is_overlapping(.{ .start_one = 3, .end_one = 4, .start_two = 55, .end_two = 57 }));
    try std.testing.expectEqual(true, hsp.is_overlapping(.{ .start_one = 3, .end_one = 4, .start_two = 20, .end_two = 28 }));
    try std.testing.expectEqual(true, hsp.is_overlapping(.{ .start_one = 3, .end_one = 4, .start_two = 20, .end_two = 25 }));
    try std.testing.expectEqual(false, hsp.is_overlapping(.{ .start_one = 3, .end_one = 4, .start_two = 20, .end_two = 24 }));
    try std.testing.expectEqual(true, hsp.is_overlapping(.{ .start_one = 0, .end_one = 100, .start_two = 0, .end_two = 100 }));
}

test "distance" {
    {
        const hsp1 = HSP{ .start_one = 1, .end_one = 1, .start_two = 25, .end_two = 25 };
        const hsp2 = HSP{ .start_one = 2, .end_one = 2, .start_two = 26, .end_two = 26 };
        try std.testing.expectEqual(@as(usize, 0), hsp1.distance_to(hsp2));
    }

    {
        const hsp1 = HSP{ .start_one = 1, .end_one = 2, .start_two = 25, .end_two = 26 };
        const hsp2 = HSP{ .start_one = 5, .end_one = 10, .start_two = 30, .end_two = 35 };
        // x: 2 -> 5 = 2 empty cells
        // y: 26 -> 30 = 3 empty cells
        // => sqrt(2*2 + 3*3)
        try std.testing.expectEqual(@as(usize, 3), hsp1.distance_to(hsp2));
        try std.testing.expectEqual(@as(usize, 3), hsp2.distance_to(hsp1));
    }
}

// test "highscore order" {
//     const allocator = std.testing.allocator;
//     var queue = std.PriorityDequeue(HSP, void, HSP.lessThanScore).init(allocator, {});
//     defer queue.deinit();

//     try queue.add(.{ .start_one = 1, .end_one = 1, .start_two = 1, .end_two = 1, .score = 55 });
//     try queue.add(.{ .start_one = 1, .end_one = 1, .start_two = 1, .end_two = 1, .score = 11 });
//     try queue.add(.{ .start_one = 1, .end_one = 1, .start_two = 1, .end_two = 1, .score = 155 });
//     try queue.add(.{ .start_one = 1, .end_one = 1, .start_two = 1, .end_two = 1, .score = -2 });

//     try std.testing.expectEqual(@as(i32, 155), queue.removeMax().score);
//     try std.testing.expectEqual(@as(i32, 55), queue.removeMax().score);
//     try std.testing.expectEqual(@as(i32, 11), queue.removeMax().score);
//     try std.testing.expectEqual(@as(i32, -2), queue.removeMax().score);
// }
