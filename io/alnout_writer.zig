const std = @import("std");

const alphabet = @import("../bio/bio.zig").alphabet;
const Sequence = @import("../sequence.zig").Sequence;
const SequenceList = @import("../sequence.zig").SequenceList;
const utils = @import("../utils.zig");
const Cigar = @import("../cigar.zig").Cigar;

pub fn AlnoutWriter(comptime A: type) type {
    return struct {
        pub fn write(writer: anytype, hits: HitList(A)) !void {
            _ = writer;
            _ = hits;

            for (hits.list.items) |hit| {
                try std.fmt.format(writer, "Query > {s}\n", .{hit.query.identifier});
            }
        }

        pub fn writeFile(path: []const u8, hits: HitList(A)) !void {
            const dir: std.fs.Dir = std.fs.cwd();
            const file: std.fs.File = try dir.createFile(path, .{});
            defer file.close();

            return try write(file.writer(), hits);
        }
    };
}

pub fn Hit(comptime A: type) type {
    return struct {
        const Self = @This();

        query: Sequence(A) = undefined,
        target: Sequence(A) = undefined,
        cigar: Cigar,

        pub fn init(query: Sequence(A), target: Sequence(A), cigar: Cigar) !Self {
            return Self{
                .query = try query.clone(),
                .target = try target.clone(),
                .cigar = try cigar.clone(),
            };
        }

        pub fn deinit(self: *Self) void {
            self.query.deinit();
            self.target.deinit();
            self.cigar.deinit();
        }
    };
}

pub fn HitList(comptime A: type) type {
    return utils.ArrayListDeinitWrapper(Hit(A));
}

test "writes correctly" {
    const DNA = alphabet.DNA;

    const allocator = std.testing.allocator;
    _ = allocator;

    var target = try Sequence(DNA).init(allocator, "DbId1", "ATCGTGAGACGATGCAAAAAATTGAGACGGATT");
    defer target.deinit();

    var query = try Sequence(DNA).init(allocator, "QryId", "GGTGAGACGTTACGCAATAAATTGAGA");
    defer query.deinit();

    var cigar = Cigar.init(allocator);
    defer cigar.deinit();

    // 2D1X8=2I1=1X4=1X9=1I
    try cigar.addWithCount(.deletion, 2);
    try cigar.addWithCount(.mismatch, 1);
    try cigar.addWithCount(.match, 8);
    try cigar.addWithCount(.insertion, 2);
    try cigar.addWithCount(.match, 1);
    try cigar.addWithCount(.mismatch, 1);
    try cigar.addWithCount(.match, 4);
    try cigar.addWithCount(.mismatch, 1);
    try cigar.addWithCount(.match, 9);
    try cigar.addWithCount(.insertion, 1);

    var hits = HitList(DNA).init(allocator);
    defer hits.deinit();

    var buffer: [1000]u8 = undefined;
    var fbs = std.io.fixedBufferStream(&buffer);
    try hits.list.append(try Hit(DNA).init(query, target, cigar));

    try AlnoutWriter(DNA).write(fbs.writer(), hits);

    var it = std.mem.split(u8, fbs.getWritten(), "\n");
    try std.testing.expectEqualStrings("Query > QryId", it.next().?);

    // try std.testing.expectEqualStrings("Query > QryId\n", fbs.getWritten());

    // Query > QryId
    //  %Id   TLen  Target
    //  81%     27   DbId1

    //  Query 28nt > QryId
    // Target 27nt > DbId1

    // Qry  1 + GGTGAGACGTTACGCAATAAATTGAGA 27
    //           ||||||||  | |||| |||||||||
    // Tgt  3 + CGTGAGACG--ATGCAAAAAATTGAGA 27

    // 27 cols, 22 ids (81.5%), 2 gaps (7.4%)
}
