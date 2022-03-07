const std = @import("std");
const print = std.debug.print;
const ArrayList = std.ArrayList;
const Allocator = std.mem.Allocator;
const builtin = @import("builtin");

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

const Args = @import("args.zig").Args;
const parseArgs = @import("args.zig").parseArgs;

const Progress = @import("progress.zig").Progress;

var progress = Progress(SearchStage){};

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
        // const NumSearchedQueriesToReport = 1;

        pub const WorkItem = struct {
            query: Sequence(DatabaseType.Alphabet),
        };

        pub const WorkerContext = struct {
            allocator: std.mem.Allocator,
            queue: *std.atomic.Queue(WorkItem),
            results: *std.atomic.Queue(Result(DatabaseType.Alphabet)),
            database: *DatabaseType,
            search_options: SearchOptions,
            search_count: *std.atomic.Atomic(usize),
            hit_count: *std.atomic.Atomic(usize),
        };

        fn entryPoint(context: *WorkerContext) !void {
            const allocator = context.allocator;
            const database = context.database;

            var search = try Search(DatabaseType).init(allocator, database, context.search_options);
            defer search.deinit();

            while (context.queue.get()) |node| {
                const query = node.data.query;

                var hits = SearchHitList(DatabaseType.Alphabet).init(allocator);
                defer hits.deinit();

                try search.search(query, &hits);

                if (hits.list.items.len > 0) {
                    const out = try context.allocator.create(std.atomic.Queue(Result(DatabaseType.Alphabet)).Node);
                    out.data = try Result(DatabaseType.Alphabet).init(allocator, query, hits);
                    context.results.put(out);

                    _ = context.hit_count.fetchAdd(hits.list.items.len, .Monotonic);
                }

                allocator.destroy(node);

                _ = context.search_count.fetchAdd(1, .Monotonic);
            }
        }
    };
}

pub fn Writer(comptime A: type) type {
    return struct {
        pub const WriterContext = struct {
            allocator: std.mem.Allocator,
            results: *std.atomic.Queue(Result(A)),
            search_finished: *std.atomic.Atomic(bool),
            written_count: *std.atomic.Atomic(usize),
            path: []const u8,
        };

        fn entryPoint(context: *WriterContext) !void {
            const dir: std.fs.Dir = std.fs.cwd();
            const file: std.fs.File = try dir.createFile(context.path, .{});
            defer file.close();

            while (true) {
                const node = context.results.get();
                if (node == null) {
                    var search_finished = context.search_finished.load(.SeqCst);
                    if (search_finished) {
                        break;
                    }

                    std.time.sleep(5 * std.time.ns_per_ms);
                    continue;
                }

                var result = node.?.data;
                defer result.deinit();

                try AlnoutWriter(A).write(file.writer(), result.query, result.hits);

                context.allocator.destroy(node.?);
                _ = context.written_count.fetchAdd(1, .Monotonic);
            }
        }
    };
}

const SearchStage = enum {
    read_database,
    analyze_database,
    index_database,
    read_queries,
    search_database,
    write_hits,
};

pub fn SearchExec(comptime A: type) type {
    const KmerLength = switch (A) {
        alphabet.DNA => 8,
        alphabet.Protein => 5,
        else => unreachable,
    };

    const databaseType = Database(A, KmerLength);
    const workerType = Worker(databaseType);
    const writerType = Writer(A);

    return struct {
        pub fn run(allocator: std.mem.Allocator, args: Args) !void {
            progress.add(.read_database, "Read database", .bytes);
            progress.add(.analyze_database, "Analyze database", .counts);
            progress.add(.index_database, "Index database", .counts);
            progress.add(.read_queries, "Read queries", .bytes);
            progress.add(.search_database, "Search database", .counts);
            progress.add(.write_hits, "Write hits", .counts);

            try progress.start();

            const filePathToDatabase = std.mem.sliceTo(&args.db, 0);
            const filePathToQuery = std.mem.sliceTo(&args.query, 0);
            const filePathToOutput = std.mem.sliceTo(&args.out, 0);
            const dir: std.fs.Dir = std.fs.cwd();

            // Read DB
            var bench_start = std.time.milliTimestamp();

            var sequences = SequenceList(A).init(allocator);
            defer sequences.deinit();

            const db_file: std.fs.File = try dir.openFile(filePathToDatabase, .{ .read = true });
            defer db_file.close();
            var db_source = std.io.StreamSource{ .file = db_file };

            var db_reader = FastaReader(A).init(allocator, &db_source);
            defer db_reader.deinit();

            progress.activate(.read_database);
            while (try db_reader.next()) |sequence| {
                try sequences.list.append(sequence);
                progress.set(.read_database, db_reader.bytes_read, db_reader.bytes_total);
            }

            // Index
            const callback = struct {
                pub fn callback(typ: databaseType.ProgressType, a: usize, b: usize) void {
                    switch (typ) {
                        .analyze => {
                            progress.activate(.analyze_database);
                            progress.set(.analyze_database, a, b);
                        },
                        .index => {
                            progress.activate(.index_database);
                            progress.set(.index_database, a, b);
                        },
                    }
                }
            }.callback;

            bench_start = std.time.milliTimestamp();
            var db = try databaseType.init(allocator, sequences.list.toOwnedSlice(), callback);
            defer db.deinit();

            // Read Query
            var queries = SequenceList(A).init(allocator);
            defer queries.deinit();

            const query_file: std.fs.File = try dir.openFile(filePathToQuery, .{ .read = true });
            defer query_file.close();
            var query_source = std.io.StreamSource{ .file = query_file };

            var query_reader = FastaReader(A).init(allocator, &query_source);
            defer query_reader.deinit();

            progress.activate(.read_queries);
            while (try query_reader.next()) |sequence| {
                try queries.list.append(sequence);
                progress.set(.read_queries, query_reader.bytes_read, query_reader.bytes_total);
            }

            // Fill queue
            var queue = std.atomic.Queue(workerType.WorkItem).init();
            var results = std.atomic.Queue(Result(A)).init();
            var search_finished = std.atomic.Atomic(bool).init(false);

            var search_count = std.atomic.Atomic(usize).init(0);
            var hit_count = std.atomic.Atomic(usize).init(0);
            var written_count = std.atomic.Atomic(usize).init(0);

            var context: workerType.WorkerContext = .{
                .allocator = allocator,
                .queue = &queue,
                .database = &db,
                .results = &results,
                .search_count = &search_count,
                .hit_count = &hit_count,
                .search_options = .{
                    .max_accepts = args.max_hits,
                    .max_rejects = args.max_rejects,
                    .min_identity = args.min_identity,
                    .strand = args.strand,
                },
            };

            var writer_context: writerType.WriterContext = .{
                .allocator = allocator,
                .results = &results,
                .search_finished = &search_finished,
                .path = filePathToOutput,
                .written_count = &written_count,
            };

            for (queries.list.items) |query| {
                const node = try context.allocator.create(std.atomic.Queue(workerType.WorkItem).Node);
                node.data = workerType.WorkItem{
                    .query = query,
                };
                queue.put(node);
            }

            // Search in threads
            bench_start = std.time.milliTimestamp();

            const worker_count = std.math.max(1, std.Thread.getCpuCount() catch 1);
            const threads = try allocator.alloc(std.Thread, worker_count);
            defer allocator.free(threads);

            var writer_thread = try std.Thread.spawn(.{}, writerType.entryPoint, .{&writer_context});

            for (threads) |*thread| {
                thread.* = try std.Thread.spawn(.{}, workerType.entryPoint, .{&context});
            }

            // wait until queue is empty
            progress.activate(.search_database);
            while (!queue.isEmpty()) {
                progress.set(.search_database, search_count.load(.SeqCst), queries.list.items.len);
                std.time.sleep(50 * std.time.ns_per_ms);
            }
            progress.set(.search_database, queries.list.items.len, queries.list.items.len); // 100%

            // wait for searches to finish
            for (threads) |thread| {
                thread.join();
            }

            // tell writer all searches are done, write remaining results and then be happy
            _ = search_finished.swap(true, .SeqCst);

            // wait until results are written
            progress.activate(.write_hits);
            while (!results.isEmpty()) {
                progress.set(.write_hits, written_count.load(.SeqCst), hit_count.load(.SeqCst));
                std.time.sleep(50 * std.time.ns_per_ms);
            }
            progress.set(.write_hits, hit_count.load(.SeqCst), hit_count.load(.SeqCst));
            writer_thread.join();

            progress.finish();
        }
    };
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();

    const allocator = switch (builtin.mode) {
        .Debug => gpa.allocator(),
        else => std.heap.c_allocator,
    };

    const args = parseArgs(allocator) catch std.os.exit(1);

    switch (args.mode) {
        .dna => try SearchExec(alphabet.DNA).run(allocator, args),
        .protein => try SearchExec(alphabet.Protein).run(allocator, args),
    }
}

test "specs" {
    std.testing.refAllDecls(@This());
}
