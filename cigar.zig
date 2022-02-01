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

    const MaxEntries = 128;

    len: usize = 0,
    buffer: [MaxEntries]Entry = undefined,

    pub fn isEmpty(self: *Self) bool {
        return self.len == 0;
    }

    pub fn clear(self: *Self) void {
        self.len = 0;
    }

    pub fn addWithCount(self: *Self, op: CigarOp, count: usize) void {
        var last_entry = if (self.len == 0) null else &self.buffer[self.len - 1];

        if (last_entry == null or last_entry.?.op != op) {
            self.buffer[self.len] = .{ .op = op, .count = count };
            self.len += 1;

            std.debug.assert(self.len < MaxEntries);
        } else {
            last_entry.?.count += count;
        }
    }

    pub fn add(self: *Self, op: CigarOp) void {
        self.addWithCount(op, 1);
    }

    pub fn appendOther(self: *Self, other: Self) void {
        for (other.entries()) |other_entry| {
            self.addWithCount(other_entry.op, other_entry.count);
        }
    }

    pub fn reverse(self: *Self) void {
        std.mem.reverse(Entry, self.buffer[0..self.len]);
    }

    pub fn toStringAlloc(self: *Self, allocator: std.mem.Allocator) ![]u8 {
        var str = std.ArrayList(u8).init(allocator);
        var buf: [128]u8 = undefined;

        for (self.entries()) |entry| {
            const xyz = try std.fmt.bufPrint(buf[0..], "{}{c}", .{ entry.count, @enumToInt(entry.op) });
            try str.appendSlice(xyz);
        }

        return str.toOwnedSlice();
    }

    pub fn identity(self: *Self) f32 {
        var letters: usize = 0;
        var matches: usize = 0;

        for (self.entries()) |entry, index| {
            const is_gap = (entry.op == .insertion or entry.op == .deletion);

            // Don't count terminal gaps towards identity calculation
            if (index == 0 and is_gap)
                continue;

            if (index == (self.len - 1) and is_gap)
                continue;

            letters += entry.count;
            matches += if (entry.op == .match) entry.count else 0;
        }

        if (letters == 0)
            return 0.0;

        return @intToFloat(f32, matches) / @intToFloat(f32, letters);
    }

    pub fn entries(self: Self) []const Entry {
        return self.buffer[0..self.len];
    }
};

test "basic add" {
    const allocator = std.testing.allocator;

    var cigar = Cigar{};

    cigar.add(CigarOp.match);
    cigar.add(CigarOp.match);
    cigar.add(CigarOp.mismatch);
    cigar.add(CigarOp.deletion);
    cigar.addWithCount(CigarOp.deletion, 2);
    cigar.add(CigarOp.insertion);
    cigar.add(CigarOp.insertion);
    cigar.add(CigarOp.match);

    const cigar_str = try cigar.toStringAlloc(allocator);
    defer allocator.free(cigar_str);

    try std.testing.expectEqualStrings("2=1X3D2I1=", cigar_str);
}

test "append other cigar" {
    const allocator = std.testing.allocator;

    var cigar = Cigar{};
    var other_cigar = Cigar{};

    // tail does not match, simple append
    {
        cigar.clear();
        cigar.add(CigarOp.match);
        cigar.add(CigarOp.match);
        cigar.add(CigarOp.mismatch);

        other_cigar.clear();
        other_cigar.add(CigarOp.match);
        other_cigar.add(CigarOp.deletion);

        cigar.appendOther(other_cigar);

        const cigar_str = try cigar.toStringAlloc(allocator);
        defer allocator.free(cigar_str);
        try std.testing.expectEqualStrings("2=1X1=1D", cigar_str);
    }

    // tail matches, expand
    {
        cigar.clear();
        cigar.add(CigarOp.match);
        cigar.add(CigarOp.match);
        cigar.add(CigarOp.mismatch);

        other_cigar.clear();
        other_cigar.add(CigarOp.mismatch);
        other_cigar.add(CigarOp.mismatch);
        other_cigar.add(CigarOp.deletion);

        cigar.appendOther(other_cigar);

        const cigar_str = try cigar.toStringAlloc(allocator);
        defer allocator.free(cigar_str);
        try std.testing.expectEqualStrings("2=3X1D", cigar_str);
    }
}
