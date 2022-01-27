const std = @import("std");

pub const Highscores = struct {
    const Self = @This();
    pub const Entry = struct {
        id: usize,
        score: isize,
    };

    allocator: std.mem.Allocator,
    entries: []Entry,

    pub fn init(allocator: std.mem.Allocator, num_entries: usize) !Self {
        std.debug.assert(num_entries > 0);

        const entries = try allocator.alloc(Entry, num_entries);
        std.mem.set(Entry, entries, Entry{ .id = 0, .score = -1 });

        return Self{
            .allocator = allocator,
            .entries = entries,
        };
    }

    pub fn add(self: *Self, id: usize, score: usize) void {
        var lowest_entry = &self.entries[self.entries.len - 1];
        if (score < lowest_entry.*.score)
            return;

        var slot: usize = for (self.entries) |entry, index| {
            if (entry.id == id) {
                // ensure score is not decreasing for existing entry
                // otherwise sorting mechanism below fails
                std.debug.assert(entry.score <= @intCast(isize, score));

                break index;
            }
        } else self.entries.len - 1;

        var entry = &self.entries[slot];
        entry.id = id;
        entry.score = @intCast(isize, score);

        // check if neighbouring entry is lower, if so swap
        while (slot >= 1) : (slot -= 1) {
            var right_entry = &self.entries[slot];
            var left_entry = &self.entries[slot - 1];

            if (left_entry.score > right_entry.score)
                break;

            swap(right_entry, left_entry);
        }
    }

    pub fn top_to_bottom(self: *Self) []const Entry {
        // skip empty entries
        var max_slot: usize = for (self.entries) |entry, index| {
            if (entry.score < 0) break index;
        } else self.entries.len; 

        return self.entries[0..max_slot];
    }

    pub fn deinit(self: *Self) void {
        self.allocator.free(self.entries);
    }

    fn swap(a: anytype, b: anytype) void {
        const tmp = a.*;
        a.* = b.*;
        b.* = tmp;
    }
};

test "test" {
    const allocator = std.testing.allocator;

    var hs = try Highscores.init(allocator, 3);
    defer hs.deinit();

    try std.testing.expectEqual(@as(usize, 0), hs.top_to_bottom().len);

    const IdA = 3;
    const IdB = 51;
    const IdC = 77;
    const IdD = 101;

    hs.add(IdA, 100);
    try std.testing.expectEqual(@as(usize, 1), hs.top_to_bottom().len);

    hs.add(IdB, 20);
    try std.testing.expectEqual(@as(usize, 2), hs.top_to_bottom().len);

    var cmp1 = [_]Highscores.Entry{ .{ .id = IdA, .score = 100 }, .{ .id = IdB, .score = 20 } };
    try std.testing.expectEqualSlices(Highscores.Entry, &cmp1, hs.top_to_bottom());

    hs.add(IdB, 150);
    try std.testing.expectEqual(@as(usize, 2), hs.top_to_bottom().len);

    hs.add(IdC, 50);
    try std.testing.expectEqual(@as(usize, 3), hs.top_to_bottom().len);

    hs.add(IdD, 30);
    try std.testing.expectEqual(@as(usize, 3), hs.top_to_bottom().len);

    var cmp2 = [_]Highscores.Entry{ .{ .id = IdB, .score = 150 }, .{ .id = IdA, .score = 100 }, .{ .id = IdC, .score = 50 } }; 
    try std.testing.expectEqualSlices(Highscores.Entry, &cmp2, hs.top_to_bottom());

    hs.add(IdC, 200);
    var cmp3 = [_]Highscores.Entry{ .{ .id = IdC, .score = 200 }, .{ .id = IdB, .score = 150 }, .{ .id = IdA, .score = 100 } };
    try std.testing.expectEqualSlices(Highscores.Entry, &cmp3, hs.top_to_bottom());
}
