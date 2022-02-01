const std = @import("std");

pub fn dup(comptime T: type, allocator: std.mem.Allocator, slice: []const T) ![]T {
    const result = try allocator.alloc(T, slice.len);
    std.mem.copy(u8, result, slice);
    return result;
}

pub fn ArrayListWithDeinit(comptime T: type) type {
    return struct {
        const Self = @This();

        allocator: std.mem.Allocator,
        items: std.ArrayList(T),

        pub fn init(allocator: std.mem.Allocator) Self {
            return Self{
                .allocator = allocator,
                .items = std.ArrayList(T).init(allocator),
            };
        }

        pub fn deinit(self: *Self) void {
            for (self.items.items) |*item|
                item.deinit();

            self.items.deinit();
        }

        pub fn append(self: *Self, item: T) !void {
            try self.items.append(item);
        }

        pub fn toOwnedSlice(self: *Self) []T {
            return self.items.toOwnedSlice();
        }

        pub fn clear(self: *Self) void {
            return self.items.clearRetainingCapacity();
        }
    };
}
