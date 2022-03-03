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
const SearchOptions = @import("search.zig").SearchOptions;
const SearchHitList = @import("search.zig").SearchHitList;

const Sequence = @import("sequence.zig").Sequence;
const SequenceList = @import("sequence.zig").SequenceList;

const FastaReader = @import("io/io.zig").FastaReader;
const AlnoutWriter = @import("io/io.zig").AlnoutWriter;

const Database = @import("database.zig").Database;
const Highscores = @import("highscores.zig").Highscores;

const bio = @import("bio/bio.zig");
const alphabet = bio.alphabet;

const CreateWorker = @import("worker.zig").CreateWorker;

const getArgs = @import("args.zig").getArgs;

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

// pub fn Worker(comptime DatabaseType: type) type {
//     return struct {
//         pub const WorkItem = struct {
//             query: Sequence(DatabaseType.Alphabet),
//         };

//         pub const WorkerContext = struct {
//             allocator: std.mem.Allocator,
//             queue: *std.atomic.Queue(WorkItem),
//             results: *std.atomic.Queue(Result(DatabaseType.Alphabet)),
//             database: *DatabaseType,
//             search_options: SearchOptions,
//         };

//         fn entryPoint(context: *WorkerContext) !void {
//             const allocator = context.allocator;
//             const database = context.database;

//             var search = try Search(DatabaseType).init(allocator, database, context.search_options);
//             defer search.deinit();

//             while (context.queue.get()) |node| {
//                 const query = node.data.query;

//                 var hits = SearchHitList(DatabaseType.Alphabet).init(allocator);
//                 defer hits.deinit();

//                 try search.search(query, &hits);

//                 if (hits.list.items.len > 0) {
//                     const out = try context.allocator.create(std.atomic.Queue(Result(DatabaseType.Alphabet)).Node);
//                     out.* = .{
//                         .data = try Result(alphabet.DNA).init(allocator, query, hits),
//                     };
//                     context.results.put(out);
//                 }

//                 allocator.destroy(node);
//             }
//         }
//     };
// }

// pub fn Writer(comptime A: type) type {
//     return struct {
//         pub const WriterContext = struct {
//             allocator: std.mem.Allocator,
//             results: *std.atomic.Queue(Result(A)),
//             search_finished: *std.atomic.Atomic(bool),
//             path: []const u8,
//         };

//         fn entryPoint(context: *WriterContext) !void {
//             const dir: std.fs.Dir = std.fs.cwd();
//             const file: std.fs.File = try dir.createFile(context.path, .{});
//             defer file.close();

//             while (true) {
//                 const node = context.results.get();
//                 if (node == null) {
//                     var search_finished = context.search_finished.load(.SeqCst);
//                     if (search_finished) {
//                         break;
//                     }

//                     std.time.sleep(5 * std.time.ns_per_ms);
//                     continue;
//                 }

//                 var result = node.?.data;
//                 defer result.deinit();

//                 try AlnoutWriter(alphabet.DNA).write(file.writer(), result.query, result.hits);

//                 context.allocator.destroy(node.?);
//             }
//         }
//     };
// }

pub fn SequenceReader(comptime A: type) type {
    return struct {
        const Self = @This();

        pub const Worker = CreateWorker(void, Sequence(A), 20);

        path: []const u8,

        pub fn handle(self: *Self, worker: *Worker) !void {
            const dir: std.fs.Dir = std.fs.cwd();
            const file: std.fs.File = try dir.openFile(self.path, .{ .read = true });
            defer file.close();

            var source = std.io.StreamSource{ .file = file };

            var reader = FastaReader(A).init(worker.allocator, &source);
            defer reader.deinit();

            while (try reader.next()) |sequence| {
                const node = try worker.allocator.create(Worker.OutputNode);
                node.data = sequence;
                worker.push(node);
            }
        }
    };
}

pub fn Indexer(comptime DatabaseType: type) type {
    return struct {
        const Self = @This();

        pub const Worker = CreateWorker(Sequence(DatabaseType.Alphabet), void, 20);

        database: *DatabaseType,

        pub fn handle(self: *Self, worker: *Worker) !void {
            var sequences = SequenceList(DatabaseType.Alphabet).init(worker.allocator);
            defer sequences.deinit();

            while (worker.pop()) |node| {
                try sequences.list.append(node.data);
                worker.allocator.destroy(node);
            }

            self.database.* = try DatabaseType.init(worker.allocator, sequences.list.toOwnedSlice());

            std.debug.print("DB ENDE \n", .{});
        }
    };
}

pub fn Searcher(comptime DatabaseType: type) type {
    return struct {
        const Self = @This();

        pub const Worker = CreateWorker(Sequence(DatabaseType.Alphabet), Result(DatabaseType.Alphabet), 20);

        database: *DatabaseType,
        search_options: SearchOptions,

        pub fn handle(self: *Self, worker: *Worker) !void {
            var search = try Search(DatabaseType).init(worker.allocator, self.database, self.search_options);
            defer search.deinit();

            while (worker.pop()) |node| {
                var hits = SearchHitList(DatabaseType.Alphabet).init(worker.allocator);
                defer hits.deinit();

                var query = node.data;
                defer query.deinit();

                // now search!
                try search.search(query, &hits);

                if (hits.list.items.len > 0) {
                    const out = try worker.allocator.create(Worker.OutputNode);
                    out.data = try Result(DatabaseType.Alphabet).init(worker.allocator, query, hits);
                    worker.push(out);
                }

                worker.allocator.destroy(node);
            }
        }
    };
}

pub fn Writer(comptime A: type) type {
    return struct {
        const Self = @This();

        pub const Worker = CreateWorker(Result(A), void, 20);

        path: []const u8,

        pub fn handle(self: *Self, worker: *Worker) !void {
            const dir: std.fs.Dir = std.fs.cwd();
            const file: std.fs.File = try dir.createFile(self.path, .{});
            defer file.close();

            while (worker.pop()) |node| {
                var result = node.data;
                defer result.deinit();

                try AlnoutWriter(A).write(file.writer(), result.query, result.hits);

                worker.allocator.destroy(node);
            }
        }
    };
}

pub fn main() !void {
    // var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    // defer _ = gpa.deinit();
    // const allocator = gpa.allocator();
    const allocator = std.heap.c_allocator;

    const A = alphabet.DNA;
    const DatabaseType = Database(A, 8);

    const args = getArgs(allocator) catch std.os.exit(1);

    const filePathToDatabase = std.mem.sliceTo(&args.db, 0);
    const filePathToQuery = std.mem.sliceTo(&args.query, 0);
    const filePathToOutput = std.mem.sliceTo(&args.out, 0);

    var database_sequences = std.atomic.Queue(Sequence(A)).init();
    defer {
        while (database_sequences.get()) |node| {
            node.data.deinit();
            allocator.destroy(node);
        }
    }

    var query_sequences = std.atomic.Queue(Sequence(A)).init();
    defer {
        while (query_sequences.get()) |node| {
            node.data.deinit();
            allocator.destroy(node);
        }
    }

    // DB
    var database_reader_worker = SequenceReader(A).Worker.init(allocator);
    defer database_reader_worker.deinit();
    database_reader_worker.output_queue = &database_sequences;

    var database_reader: SequenceReader(A) = .{
        .path = filePathToDatabase,
    };
    try database_reader_worker.run(SequenceReader(A), &database_reader);

    // Query
    var query_reader_worker = SequenceReader(A).Worker.init(allocator);
    defer query_reader_worker.deinit();
    query_reader_worker.output_queue = &query_sequences;

    var query_reader: SequenceReader(A) = .{
        .path = filePathToQuery,
    };
    try query_reader_worker.run(SequenceReader(A), &query_reader);

    // DB indexer
    var database: DatabaseType = undefined;
    defer database.deinit();

    var indexer_worker = Indexer(DatabaseType).Worker.init(allocator);
    defer indexer_worker.deinit();
    indexer_worker.input_queue = &database_sequences;

    try indexer_worker.addDependency(.{
        .state = &database_reader_worker.state,
        .kind = .upstream,
    });

    var indexer: Indexer(DatabaseType) = .{
        .database = &database,
    };
    try indexer_worker.run(Indexer(DatabaseType), &indexer);

    // Searching
    var results = std.atomic.Queue(Result(A)).init();
    defer {
        while (results.get()) |node| {
            node.data.deinit();
            allocator.destroy(node);
        }
    }

    const worker_count = std.math.max(1, std.Thread.getCpuCount() catch 1);

    var searcher_workers = try allocator.alloc(Searcher(DatabaseType).Worker, worker_count);
    defer {
        for (searcher_workers) |*searcher_worker|
            searcher_worker.deinit();
        allocator.free(searcher_workers);
    }

    var searchers = try allocator.alloc(Searcher(DatabaseType), worker_count);
    defer allocator.free(searchers);

    for (searcher_workers) |*searcher_worker, index| {
        searcher_worker.* = Searcher(DatabaseType).Worker.init(allocator);
        searcher_worker.input_queue = &query_sequences;
        searcher_worker.output_queue = &results;

        try searcher_worker.addDependency(.{
            .state = &indexer_worker.state,
            .kind = .upstream,
        });

        searchers[index] = .{
            .database = &database,
            .search_options = .{},
        };

        try searcher_worker.run(Searcher(DatabaseType), &searchers[index]);
    }

    // Writer
    var writer_worker = Writer(A).Worker.init(allocator);
    defer writer_worker.deinit();
    writer_worker.input_queue = &results;

    for (searcher_workers) |*searcher_worker| {
        try writer_worker.addDependency(.{
            .state = &searcher_worker.state,
            .kind = .consumer,
        });
    }

    var writer: Writer(A) = .{
        .path = filePathToOutput,
    };
    try writer_worker.run(Writer(A), &writer);

    writer_worker.join();

    // const databaseType = Database(alphabet.DNA, 8);
    // const workerType = Worker(databaseType);
    // const writerType = Writer(alphabet.DNA);

    // var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    // defer _ = gpa.deinit();
    // const allocator = gpa.allocator();
    // // const allocator = std.heap.c_allocator;

    // const args = getArgs(allocator) catch std.os.exit(1);

    // // Read DB
    // var bench_start = std.time.milliTimestamp();

    // var sequences = SequenceList(alphabet.DNA).init(allocator);
    // defer sequences.deinit();

    // const db_file: std.fs.File = try dir.openFile(filePathToDatabase, .{ .read = true });
    // defer db_file.close();
    // var db_source = std.io.StreamSource{ .file = db_file };

    // var db_reader = FastaReader(alphabet.DNA).init(allocator, &db_source);
    // defer db_reader.deinit();

    // while (try db_reader.next()) |sequence| {
    //     try sequences.list.append(sequence);
    // }

    // // Read Query
    // var queries = SequenceList(alphabet.DNA).init(allocator);
    // defer queries.deinit();

    // const query_file: std.fs.File = try dir.openFile(filePathToQuery, .{ .read = true });
    // defer query_file.close();
    // var query_source = std.io.StreamSource{ .file = query_file };

    // var query_reader = FastaReader(alphabet.DNA).init(allocator, &query_source);
    // defer query_reader.deinit();

    // while (try query_reader.next()) |sequence| {
    //     try queries.list.append(sequence);
    // }

    // print("Reading took {}ms\n", .{std.time.milliTimestamp() - bench_start});

    // bench_start = std.time.milliTimestamp();
    // var db = try databaseType.init(allocator, sequences.list.toOwnedSlice());
    // defer db.deinit();

    // print("Indexing took {}ms\n", .{std.time.milliTimestamp() - bench_start});

    // // Fill queue
    // var queue = std.atomic.Queue(workerType.WorkItem).init();
    // var results = std.atomic.Queue(Result(alphabet.DNA)).init();
    // var search_finished = std.atomic.Atomic(bool).init(false);

    // var context: workerType.WorkerContext = .{ .allocator = allocator, .queue = &queue, .database = &db, .results = &results, .search_options = .{
    //     .max_accepts = args.max_hits,
    //     .max_rejects = args.max_rejects,
    //     .min_identity = args.min_identity,
    // } };

    // var writer_context: writerType.WriterContext = .{
    //     .allocator = allocator,
    //     .results = &results,
    //     .search_finished = &search_finished,
    //     .path = filePathToOutput,
    // };

    // for (queries.list.items) |query| {
    //     const node = try context.allocator.create(std.atomic.Queue(workerType.WorkItem).Node);
    //     node.* = .{
    //         .data = workerType.WorkItem{
    //             .query = query,
    //         },
    //     };
    //     queue.put(node);
    // }

    // // Search in threads
    // bench_start = std.time.milliTimestamp();

    // const worker_count = std.math.max(1, std.Thread.getCpuCount() catch 1);
    // const threads = try allocator.alloc(std.Thread, worker_count);
    // defer allocator.free(threads);

    // var writer_thread = try std.Thread.spawn(.{}, writerType.entryPoint, .{&writer_context});

    // for (threads) |*thread| {
    //     thread.* = try std.Thread.spawn(.{}, workerType.entryPoint, .{&context});
    // }

    // // wait for search to finish
    // for (threads) |thread| {
    //     thread.join();
    // }

    // // wait for write to finish
    // _ = search_finished.swap(true, .SeqCst);
    // writer_thread.join();

    // print("Searching took {}ms\n", .{std.time.milliTimestamp() - bench_start});
}
