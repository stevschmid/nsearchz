const std = @import("std");
const print = std.debug.print;
const ArrayList = std.ArrayList;
const Allocator = std.mem.Allocator;

const utils = @import("utils.zig");
const ea = @import("extend_align.zig");
const ba = @import("banded_align.zig");

const Cigar = @import("cigar.zig").Cigar;
const CigarOp = @import("cigar.zig").CigarOp;

const Search = @import("search.zig").Search;
const SearchHitList = @import("search.zig").SearchHitList;

const Sequence = @import("sequence.zig").Sequence;
const SequenceList = @import("sequence.zig").SequenceList;

const FastaReader = @import("fasta_reader.zig").FastaReader;
const Database = @import("database.zig").Database;
const Highscores = @import("highscores.zig").Highscores;

const bio = @import("bio/bio.zig");
const alphabet = bio.alphabet;

pub fn Worker(comptime DatabaseType: type) type {
    return struct {
        pub const WorkItem = struct {
            query: Sequence(DatabaseType.Alphabet),
        };

        pub const WorkerContext = struct {
            allocator: std.mem.Allocator,
            queue: *std.atomic.Queue(WorkItem),
            database: *DatabaseType,
        };

        fn entryPoint(context: *WorkerContext) !void {
            const allocator = context.allocator;
            const database = context.database;

            var search = try Search(DatabaseType).init(allocator, database, .{});
            defer search.deinit();
            _ = search;

            while (context.queue.get()) |node| {
                const query = node.data.query;

                var hits = SearchHitList.init(allocator);
                defer hits.deinit();
                try search.search(query, &hits);

                for (hits.list.items) |hit| {
                    std.debug.print("Hello {s}\n", .{hit.cigar.str()});
                }

                allocator.destroy(node);
            }
        }
    };
}

pub fn main() !void {
    const databaseType = Database(alphabet.DNA, 8);
    const workerType = Worker(databaseType);

    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    var bench_start = std.time.milliTimestamp();

    var arg_it = std.process.args();

    // skip my own exe name
    _ = arg_it.skip();

    const filePathToDatabase = (try arg_it.next(allocator) orelse {
        print("Expected first argument to be path to database file\n", .{});
        return error.InvalidArgs;
    });
    defer allocator.free(filePathToDatabase);

    const filePathToQuery = (try arg_it.next(allocator) orelse {
        print("Expected second argument to be path to query file\n", .{});
        return error.InvalidArgs;
    });
    defer allocator.free(filePathToQuery);

    // Read DB
    var sequences = SequenceList(alphabet.DNA).init(allocator);
    defer sequences.deinit();
    try FastaReader(alphabet.DNA).readFile(filePathToDatabase, &sequences);

    // Read Query
    var queries = SequenceList(alphabet.DNA).init(allocator);
    defer queries.deinit();
    try FastaReader(alphabet.DNA).readFile(filePathToQuery, &queries);

    print("Reading took {}ms\n", .{std.time.milliTimestamp() - bench_start});
    bench_start = std.time.milliTimestamp();
    var db = try databaseType.init(allocator, sequences.list.toOwnedSlice());
    defer db.deinit();
    print("Indexing took {}ms\n", .{std.time.milliTimestamp() - bench_start});

    // searching

    var queue = std.atomic.Queue(workerType.WorkItem).init();
    var context: workerType.WorkerContext = .{
        .allocator = allocator,
        .queue = &queue,
        .database = &db,
    };

    for (queries.list.items) |query| {
        const node = try context.allocator.create(std.atomic.Queue(workerType.WorkItem).Node);
        node.* = .{
            .data = workerType.WorkItem{
                .query = query,
            },
        };
        queue.put(node);
    }

    var thread = try std.Thread.spawn(.{}, workerType.entryPoint, .{&context});
    thread.join();
}
