const std = @import("std");

const HighscoreEntry = struct {
    id: usize,
    score: isize,
};

const Highscores = struct {
    const Self = @This();

    allocator: std.mem.Allocator,
    entries: []HighscoreEntry,

    pub fn init(allocator: std.mem.Allocator, num_entries: usize) !Self {
        std.debug.assert(num_entries > 0);

        const entries = try allocator.alloc(HighscoreEntry, num_entries);
        std.mem.set(HighscoreEntry, entries, HighscoreEntry{ .id = 0, .score = -1 });

        return Self{
            .allocator = allocator,
            .entries = entries,
        };
    }

    pub fn add(self: *Self, id: usize, score: usize) void {
        var lowest_entry = &self.entries[0];
        if (score < lowest_entry.*.score)
            return;

        var slot: usize = for (self.entries) |entry, index| {
            if (entry.id == id) {
                // ensure score is not decreasing for existing entry
                // otherwise sorting mechanism below fails
                std.debug.assert(entry.score <= @intCast(isize, score));

                break index;
            }
        } else 0;

        var entry = &self.entries[slot];
        entry.id = id;
        entry.score = @intCast(isize, score);

        // check if neighbouring entry is higher, if so swap
        while (slot + 1 < self.entries.len) : (slot += 1) {
            var this_entry = &self.entries[slot];
            var next_entry = &self.entries[slot + 1];

            if (next_entry.score > this_entry.score)
                break;

            swap(this_entry, next_entry);
        }
    }

    pub fn result(self: *Self) ?[]const HighscoreEntry {
        // skip empty entries
        for(self.entries) |entry, index| {
            if (entry.score >= 0) {
                return self.entries[index..];
            }
        } else return null;

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

    try std.testing.expect(hs.result() == null);

    const IdA = 3;
    const IdB = 51;
    const IdC = 77;
    const IdD = 101;

    hs.add(IdA, 100);
    try std.testing.expectEqual(@as(usize, 1), hs.result().?.len);

    hs.add(IdB, 20);
    try std.testing.expectEqual(@as(usize, 2), hs.result().?.len);

    hs.add(IdB, 150);
    try std.testing.expectEqual(@as(usize, 2), hs.result().?.len);

    hs.add(IdC, 50);
    try std.testing.expectEqual(@as(usize, 3), hs.result().?.len);

    hs.add(IdD, 30);
    try std.testing.expectEqual(@as(usize, 3), hs.result().?.len);

    var cmp = [_]HighscoreEntry{ .{ .id = IdC, .score = 50 }, .{ .id = IdA, .score = 100 }, .{ .id = IdB, .score = 150 } }; 
    try std.testing.expectEqualSlices(HighscoreEntry, &cmp, hs.result().?);

    hs.add(IdC, 200);
    cmp = [_]HighscoreEntry{ .{ .id = IdA, .score = 100 }, .{ .id = IdB, .score = 150 }, .{ .id = IdC, .score = 200 } };
    try std.testing.expectEqualSlices(HighscoreEntry, &cmp, hs.result().?);
}
