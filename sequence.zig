const std = @import("std");
const utils = @import("utils.zig");

pub fn Sequence(comptime A: type) type {
    return struct {
        const Self = @This();

        allocator: std.mem.Allocator,
        identifier: []u8,
        data: []u8,

        pub fn init(allocator: std.mem.Allocator, identifier: []const u8, data: []const u8) !Self {
            return Self{
                .allocator = allocator,
                .identifier = try utils.dup(u8, allocator, identifier),
                .data = try utils.dup(u8, allocator, data),
            };
        }

        pub fn deinit(self: *Self) void {
            self.allocator.free(self.identifier);
            self.allocator.free(self.data);
        }

        pub fn matches(self: Self, other: Self) bool {
            if (self.data.len != other.data.len) {
                return false;
            }

            return for (self.data) |_, index| {
                if (self.data[index] != other.data[index]) break false;
            } else true;
        }

        pub fn clone(self: *Self) !Self {
            return Self{
                .allocator = self.allocator,
                .identifier = try utils.dup(u8, self.allocator, self.identifier),
                .data = try utils.dup(u8, self.allocator, self.data),
            };
        }

        pub fn reverse(self: *Self) void {
            // reverse
            std.mem.reverse(u8, self.data);
        }

        pub fn complement(self: *Self) void {
            // complement
            for (self.data) |*letter| {
                letter.* = A.complement(letter.*);
            }
        }
    };
}

pub fn SequenceList(comptime A: type) type {
    return struct {
        const Self = @This();

        allocator: std.mem.Allocator,
        list: std.ArrayList(Sequence(A)),

        pub fn init(allocator: std.mem.Allocator) Self {
            return Self{
                .allocator = allocator,
                .list = std.ArrayList(Sequence(A)).init(allocator),
            };
        }

        pub fn append(self: *Self, identifier: []const u8, data: []const u8) !void {
            const seq = try Sequence(A).init(self.allocator, identifier, data);
            try self.list.append(seq);
        }

        pub fn toOwnedSlice(self: *Self) []Sequence(A) {
            return self.list.toOwnedSlice();
        }

        pub fn deinit(self: *Self) void {
            for (self.list.items) |*seq|
                seq.deinit();

            self.list.deinit();
        }
    };
}
