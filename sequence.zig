const std = @import("std");
const utils = @import("utils.zig");

pub fn Sequence(comptime A: type) type {
    return struct {
        const Self = @This();

        identifier: []const u8,
        data: []const u8,

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

pub fn SequenceStore(comptime A: type) type {
    return struct {
        const Self = @This();

        allocator: std.mem.Allocator,
        store: std.ArrayList(Sequence(A)),

        pub fn init(allocator: std.mem.Allocator) Self {
            return Self{
                .allocator = allocator,
                .store = std.ArrayList(Sequence(A)).init(allocator),
            };
        }

        pub fn appendSeq(self: *Self, seq: Sequence(A)) !void {
            return try self.append(seq.identifier, seq.data);
        }

        pub fn append(self: *Self, identifier_: []const u8, data_: []const u8) !void {
            const identifier = try utils.dup(u8, self.allocator, identifier_);
            const data = try utils.dup(u8, self.allocator, data_);

            const sequence: Sequence(A) = .{ .identifier = identifier, .data = data };

            try self.store.append(sequence);
        }

        pub fn deinit(self: *Self) void {
            for (self.store.items) |*seq| {
                self.allocator.free(seq.identifier);
                self.allocator.free(seq.data);
            }

            self.store.deinit();
        }

        pub fn sequences(self: *Self) []Sequence(A) {
            return self.store.items;
        }
    };
}
