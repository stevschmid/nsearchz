const std = @import("std");
const utils = @import("utils.zig");

pub fn Sequence(comptime A: type) type {
    return struct {
        const Self = @This();

        allocator: std.mem.Allocator,

        identifier: []const u8,
        data: []const u8,

        pub fn init(allocator: std.mem.Allocator, identifier: []const u8, data: []const u8) !Self {
            var sanitizedData = try utils.dup(u8, allocator, data);

            // convert lower case to upper case
            for (sanitizedData) |*letter| {
                if (letter.* >= 'a' and letter.* <= 'z') {
                    letter.* -= ('a' - 'A');
                }
            }


            return Self {
                .allocator = allocator,
                .identifier = try utils.dup(u8, allocator, identifier),
                .data = sanitizedData,
            };
        }

        pub fn deinit(self: *Self) void {
            self.allocator.free(self.identifier);
            self.allocator.free(self.data);
        }

        pub fn matches(self: *const Self, other: Self) bool {
            if (self.data.len != other.data.len) {
                return false;
            }

            return for (self.data) |_, index| {
                if (self.data[index] != other.data[index]) break false;
            } else true;
        }

        pub fn complementAlloc(self: *const Self) !Self {
            var compl_data = try utils.dup(u8, self.allocator, self.data);

            for (compl_data) |*letter| {
                letter.* = A.complement(letter.*);
            }

            return Sequence(A) {
                .allocator = self.allocator,
                .identifier = try utils.dup(u8, self.allocator, self.identifier),
                .data = compl_data,
            };
        }
    };
}
