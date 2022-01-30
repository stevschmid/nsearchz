const std = @import("std");
const utils = @import("utils.zig");

pub fn Sequence(comptime A: type) type {
    return struct {
        const Self = @This();

        allocator: std.mem.Allocator,

        identifier: []const u8,
        data: []const u8,

        pub fn init(allocator: std.mem.Allocator, identifier_: []const u8, data_: []const u8) !Self {
            var identifier = try utils.dup(u8, allocator, identifier_);
            var data = try utils.dup(u8, allocator, data_);

            // convert lower case to upper case
            for (data) |*letter| {
                if (letter.* >= 'a' and letter.* <= 'z') {
                    letter.* -= ('a' - 'A');
                }
            }

            return Self{
                .allocator = allocator,
                .identifier = identifier,
                .data = data,
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

        pub fn complementAlloc(allocator: std.mem.Allocator, self: *const Self) !Self {
            var compl_identifier = try utils.dup(u8, allocator, self.identifier);
            var compl_data = try utils.dup(u8, allocator, self.data);

            for (compl_data) |*letter| {
                letter.* = A.complement(letter.*);
            }

            return Sequence(A){
                .allocator = .allocator,
                .identifier = compl_identifier,
                .data = compl_data,
            };
        }
    };
}
