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

    pub fn clear(self: *Self) void {
        self.entries.clearRetainingCapacity();
    }

    pub fn add(self: *Self, op: CigarOp) !void {
        var last_entry = if (self.entries.items.len == 0) null else &self.entries.items[self.entries.items.len - 1];
        if (last_entry == null or last_entry.?.op != op) {
            try self.entries.append(Entry{ .op = op, .count = 1 });
        } else {
            last_entry.?.count += 1;
        }
    }

    pub fn reverse(self: *Self) void {
        std.mem.reverse(Entry, self.entries.items);
    }

    pub fn toStringAlloc(self: *Self, allocator: std.mem.Allocator) ![]u8 {
        var str = std.ArrayList(u8).init(allocator);
        var buf: [128]u8 = undefined;

        for (self.entries.items) |entry| {
            const xyz = try std.fmt.bufPrint(buf[0..], "{}{c}", .{ entry.count, @enumToInt(entry.op) });
            try str.appendSlice(xyz);
        }

        return str.toOwnedSlice();
    }

    pub fn deinit(self: *Self) void {
        self.entries.deinit();
    }
};

// TODO: tests
