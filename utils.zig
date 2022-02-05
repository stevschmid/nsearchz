const std = @import("std");

pub fn dup(comptime T: type, allocator: std.mem.Allocator, slice: []const T) ![]T {
    const result = try allocator.alloc(T, slice.len);
    std.mem.copy(u8, result, slice);
    return result;
}

pub fn ArrayListDeinitWrapper(comptime T: type) type {
    return struct {
        const Self = @This();

        allocator: std.mem.Allocator,
        list: std.ArrayList(T),

        pub fn init(allocator: std.mem.Allocator) Self {
            return Self{
                .allocator = allocator,
                .list = std.ArrayList(T).init(allocator),
            };
        }

        pub fn clone(self: Self) !Self {
            var cloned = Self.init(self.allocator);

            for (self.list.items) |item|
                try cloned.list.append(try item.clone());

            return cloned;
        }

        pub fn deinit(self: *Self) void {
            for (self.list.items) |*item|
                item.deinit();

            self.list.deinit();
        }
    };
}
