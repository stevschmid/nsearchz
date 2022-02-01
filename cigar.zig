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

    entries: std.ArrayList(Entry),

    pub fn init(allocator: std.mem.Allocator) Self {
        return Self{
            .entries = std.ArrayList(Entry).init(allocator),
        };
    }

    pub fn isEmpty(self: *Self) bool {
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

    pub fn appendOther(self: *Self, other: Self) !void {
        for (other.entries.items) |other_entry| {
            try self.addWithCount(other_entry.op, other_entry.count);
        }
    }

    pub fn reverse(self: *Self) void {
        std.mem.reverse(Entry, self.entries.items);
    }

    pub fn toStringAlloc(self: Self, allocator: std.mem.Allocator) ![]u8 {
        var str = std.ArrayList(u8).init(allocator);
        var buf: [128]u8 = undefined;

        for (self.entries.items) |entry| {
            const xyz = try std.fmt.bufPrint(buf[0..], "{}{c}", .{ entry.count, @enumToInt(entry.op) });
            try str.appendSlice(xyz);
        }

        return str.toOwnedSlice();
    }

    pub fn identity(self: *Self) f32 {
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

    try cigar.add(CigarOp.match);
    try cigar.add(CigarOp.match);
    try cigar.add(CigarOp.mismatch);
    try cigar.add(CigarOp.deletion);
    try cigar.addWithCount(CigarOp.deletion, 2);
    try cigar.add(CigarOp.insertion);
    try cigar.add(CigarOp.insertion);
    try cigar.add(CigarOp.match);

    const cigar_str = try cigar.toStringAlloc(allocator);
    defer allocator.free(cigar_str);

    try std.testing.expectEqualStrings("2=1X3D2I1=", cigar_str);
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
        try cigar.add(CigarOp.match);
        try cigar.add(CigarOp.match);
        try cigar.add(CigarOp.mismatch);

        other_cigar.clear();
        try other_cigar.add(CigarOp.match);
        try other_cigar.add(CigarOp.deletion);

        try cigar.appendOther(other_cigar);

        const cigar_str = try cigar.toStringAlloc(allocator);
        defer allocator.free(cigar_str);
        try std.testing.expectEqualStrings("2=1X1=1D", cigar_str);
    }

    // tail matches, expand
    {
        cigar.clear();
        try cigar.add(CigarOp.match);
        try cigar.add(CigarOp.match);
        try cigar.add(CigarOp.mismatch);

        other_cigar.clear();
        try other_cigar.add(CigarOp.mismatch);
        try other_cigar.add(CigarOp.mismatch);
        try other_cigar.add(CigarOp.deletion);

        try cigar.appendOther(other_cigar);

        const cigar_str = try cigar.toStringAlloc(allocator);
        defer allocator.free(cigar_str);
        try std.testing.expectEqualStrings("2=3X1D", cigar_str);
    }
}
