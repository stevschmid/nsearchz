const std = @import("std");

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
