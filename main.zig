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

const FastaReader = @import("io/io.zig").FastaReader;
const AlnoutWriter = @import("io/io.zig").AlnoutWriter;

const Database = @import("database.zig").Database;
const Highscores = @import("highscores.zig").Highscores;

const bio = @import("bio/bio.zig");
const alphabet = bio.alphabet;

pub fn Result(comptime A: type) type {
    return struct {
        const Self = @This();

        allocator: std.mem.Allocator,
        query: Sequence(A),
        hits: SearchHitList(A),

        pub fn init(allocator: std.mem.Allocator, query: Sequence(A), hits: SearchHitList(A)) !Self {
            return Self{
                .allocator = allocator,
                .query = try query.clone(),
                .hits = try hits.clone(),
            };
        }

        pub fn deinit(self: *Self) void {
            self.query.deinit();
            self.hits.deinit();
        }
    };
}

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

            while (context.queue.get()) |node| {
                const query = node.data.query;

                var hits = SearchHitList(DatabaseType.Alphabet).init(allocator);
                defer hits.deinit();

                try search.search(query, &hits);

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

    const filePathToOutput = (try arg_it.next(allocator) orelse {
        print("Expected third argument to be path to output file\n", .{});
        return error.InvalidArgs;
    });
    defer allocator.free(filePathToOutput);

    // Read DB
    var bench_start = std.time.milliTimestamp();

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

    // Fill queue
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

    // Search in threads
    bench_start = std.time.milliTimestamp();

    const worker_count = std.math.max(1, std.Thread.getCpuCount() catch 1);
    const threads = try allocator.alloc(std.Thread, worker_count);
    defer allocator.free(threads);

    for (threads) |*thread| {
        thread.* = try std.Thread.spawn(.{}, workerType.entryPoint, .{&context});
    }

    for (threads) |thread| {
        thread.join();
    }

    print("Searching took {}ms\n", .{std.time.milliTimestamp() - bench_start});

    //     var search = try Search(databaseType).init(allocator, &db, .{});
    //     defer search.deinit();

    //     var results = utils.ArrayListDeinitWrapper(Result(alphabet.DNA)).init(allocator);
    //     defer results.deinit();

    //     bench_start = std.time.milliTimestamp();
    //     for (queries.list.items) |query| {
    //         var hits = SearchHitList(alphabet.DNA).init(allocator);
    //         defer hits.deinit();

    //         try search.search(query, &hits);

    //         try results.list.append(try Result(alphabet.DNA).init(allocator, query, hits));
    //     }

    //     const dir: std.fs.Dir = std.fs.cwd();
    //     const file: std.fs.File = try dir.createFile(filePathToOutput, .{});
    //     defer file.close();

    //     bench_start = std.time.milliTimestamp();
    //     for (results.list.items) |result| {
    //         try AlnoutWriter(alphabet.DNA).write(file.writer(), result.query, result.hits);
    //     }
    //     print("Writing took {}ms\n", .{std.time.milliTimestamp() - bench_start});
}
