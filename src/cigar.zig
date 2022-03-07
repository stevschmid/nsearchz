const std = @import("std");

pub const CigarOp = enum(u8) {
    match = '=',
    mismatch = 'X',
    deletion = 'D',
    insertion = 'I',
};

pub const Cigar = struct {
    const Self = @This();
    const Entry = struct {
        op: CigarOp,
        count: usize = 0,
    };

    allocator: std.mem.Allocator,
    entries: std.ArrayList(Entry),

    pub fn init(allocator: std.mem.Allocator) Self {
        return Self{
            .allocator = allocator,
            .entries = std.ArrayList(Entry).init(allocator),
        };
    }

    pub fn clone(self: Self) !Self {
        var cloned = Self.init(self.allocator);
        try cloned.entries.appendSlice(self.entries.items);
        return cloned;
    }

    pub fn isEmpty(self: Self) bool {
        return self.entries.items.len == 0;
    }

    pub fn clear(self: *Self) void {
        self.entries.clearRetainingCapacity();
    }

    pub fn addWithCount(self: *Self, op: CigarOp, count: usize) !void {
        var last_entry = if (self.entries.items.len == 0) null else &self.entries.items[self.entries.items.len - 1];

        if (last_entry == null or last_entry.?.op != op) {
            try self.entries.append(Entry{ .op = op, .count = count });
        } else {
            last_entry.?.count += count;
        }
    }

    pub fn add(self: *Self, op: CigarOp) !void {
        try self.addWithCount(op, 1);
    }

    pub fn addFromString(self: *Self, s: []const u8) !void {
        var counts = std.mem.tokenize(u8, s, "=XDI");
        var ops = std.mem.tokenize(u8, s, "0123456789");

        while (counts.next()) |count| {
            var op = ops.next();
            if (op == null)
                break;

            try self.addWithCount(@intToEnum(CigarOp, op.?[0]), try std.fmt.parseInt(usize, count, 10));
        }
    }

    pub fn appendOther(self: *Self, other: Self) !void {
        for (other.entries.items) |other_entry| {
            try self.addWithCount(other_entry.op, other_entry.count);
        }
    }

    pub const Iterator = struct {
        cigar: *const Self,
        index: usize = 0,
        count: usize = 0,

        pub fn next(it: *Iterator) ?CigarOp {
            if (it.isEmpty())
                return null;

            var entry = it.cigar.entries.items[it.index];

            it.count += 1;
            if (it.count >= entry.count) {
                // next element
                it.count = 0;
                it.index += 1;
            }

            return entry.op;
        }

        pub fn isEmpty(it: *Iterator) bool {
            if (it.cigar.entries.items.len == 0)
                return true;

            if (it.index >= it.cigar.entries.items.len)
                return true;

            return false;
        }
    };

    pub fn iterator(self: *const Self) Iterator {
        return Iterator{ .cigar = self };
    }

    pub fn reverse(self: *Self) void {
        std.mem.reverse(Entry, self.entries.items);
    }

    pub fn str(self: Self) []const u8 {
        // let's try this nifty local buffer thingy by idtech idStr
        const S = struct {
            threadlocal var buffer: [10_000]u8 = undefined;
        };

        var fbs = std.io.fixedBufferStream(&S.buffer);
        var buf: [128]u8 = undefined;

        for (self.entries.items) |entry| {
            const xyz = std.fmt.bufPrint(&buf, "{}{c}", .{ entry.count, @enumToInt(entry.op) }) catch &[_]u8{};
            _ = fbs.write(xyz) catch 0;
        }

        return fbs.getWritten();
    }

    pub fn identity(self: Self) f32 {
        var letters: usize = 0;
        var matches: usize = 0;

        for (self.entries.items) |entry, index| {
            const is_gap = (entry.op == .insertion or entry.op == .deletion);

            // Don't count terminal gaps towards identity calculation
            if (index == 0 and is_gap)
                continue;

            if (index == self.entries.items.len - 1 and is_gap)
                continue;

            letters += entry.count;
            matches += if (entry.op == .match) entry.count else 0;
        }

        if (letters == 0)
            return 0.0;

        return @intToFloat(f32, matches) / @intToFloat(f32, letters);
    }

    pub fn deinit(self: *Self) void {
        self.entries.deinit();
    }
};

test "basic add" {
    const allocator = std.testing.allocator;

    var cigar = Cigar.init(allocator);
    defer cigar.deinit();

    try cigar.add(.match);
    try cigar.add(.match);
    try cigar.add(.mismatch);
    try cigar.add(.deletion);
    try cigar.addWithCount(.deletion, 2);
    try cigar.add(.insertion);
    try cigar.add(.insertion);
    try cigar.add(.match);

    try std.testing.expectEqualStrings("2=1X3D2I1=", cigar.str());
}

test "pop op" {
    const allocator = std.testing.allocator;

    var cigar = Cigar.init(allocator);
    defer cigar.deinit();

    try cigar.add(.match);
    try cigar.add(.mismatch);
    try cigar.addWithCount(.insertion, 2);
    try cigar.add(.match);

    var it = cigar.iterator();

    try std.testing.expectEqual(CigarOp.match, it.next().?);
    try std.testing.expectEqual(CigarOp.mismatch, it.next().?);
    try std.testing.expectEqual(CigarOp.insertion, it.next().?);
    try std.testing.expectEqual(CigarOp.insertion, it.next().?);
    try std.testing.expectEqual(CigarOp.match, it.next().?);
    try std.testing.expect(it.next() == null);
}

test "append other cigar" {
    const allocator = std.testing.allocator;

    var cigar = Cigar.init(allocator);
    defer cigar.deinit();

    var other_cigar = Cigar.init(allocator);
    defer other_cigar.deinit();
    // tail does not match, simple append
    {
        cigar.clear();
        try cigar.add(.match);
        try cigar.add(.match);
        try cigar.add(.mismatch);

        other_cigar.clear();
        try other_cigar.add(.match);
        try other_cigar.add(.deletion);

        try cigar.appendOther(other_cigar);

        try std.testing.expectEqualStrings("2=1X1=1D", cigar.str());
    }

    // tail matches, expand
    {
        cigar.clear();
        try cigar.add(.match);
        try cigar.add(.match);
        try cigar.add(.mismatch);

        other_cigar.clear();
        try other_cigar.add(.mismatch);
        try other_cigar.add(.mismatch);
        try other_cigar.add(.deletion);

        try cigar.appendOther(other_cigar);

        try std.testing.expectEqualStrings("2=3X1D", cigar.str());
    }
}

test "addFromString" {
    const allocator = std.testing.allocator;

    var cigar = Cigar.init(allocator);
    defer cigar.deinit();

    try cigar.addFromString("5=2I3X25=3D");
    try std.testing.expectEqualStrings("5=2I3X25=3D", cigar.str());
}
