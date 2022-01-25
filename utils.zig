const std = @import("std");

pub fn dup(comptime T: type, allocator: std.mem.Allocator, slice: []const T) ![]T {
    const result = try allocator.alloc(T, slice.len);
    std.mem.copy(u8, result, slice);
    return result;
}

