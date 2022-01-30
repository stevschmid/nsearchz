const std = @import("std");
const alphabet = @import("bio/bio.zig").alphabet;

const Database = @import("database.zig").Database;
const Sequence = @import("sequence.zig").Sequence;
const Highscores = @import("highscores.zig").Highscores;

const ea = @import("extend_align.zig");
const ba = @import("banded_align.zig");

// const Strand = enum {
//     Plus,
//     Minus,
//     Both,
// };

const SearchOptions = struct {
    max_accepts: u32 = 1,
    max_rejects: u32 = 16,
    min_identity: f32 = 0.8,
    // strand: Strand = Strand.Plus,
};

pub fn Search(comptime DatabaseType: type) type {
    return struct {
        const Self = @This();
        const A = DatabaseType.kmerInfo.Alphabet;

        const DefaultMinHspLength = 16;
        const DefaultMaxHSPJoinDistance = 16;
        const kmerInfo = DatabaseType.kmerInfo;

        allocator: std.mem.Allocator,
        options: SearchOptions,
        database: DatabaseType,
        extend_align: ea.ExtendAlign(A),
        banded_align: ba.BandedAlign(A),
        kmers: []kmerInfo.KmerType,

        pub fn init(allocator: std.mem.Allocator, database: DatabaseType, options: SearchOptions) !Self {
            return Self{
                .allocator = allocator,
                .database = database,
                .options = options,
                .extend_align = ea.ExtendAlign(A).init(allocator, .{}),
                .banded_align = ba.BandedAlign(A).init(allocator, .{}),
                .kmers = &[_]kmerInfo.KmerType{},
            };
        }

        pub fn process(self: *Self, query: Sequence(A)) !void {
            const min_hsp_length = std.math.min(DefaultMinHspLength, query.data.len / 2);
            const max_hsp_join_distance = DefaultMaxHSPJoinDistance;
            _ = min_hsp_length;
            _ = max_hsp_join_distance;

            var highscores = try Highscores.init(self.allocator, self.options.max_accepts + self.options.max_rejects);
            defer highscores.deinit();

            _ = query;
        }

        pub fn deinit(self: *Self) void {
            self.extend_align.deinit();
            self.banded_align.deinit();
        }
    };
}

test "check" {
    const allocator = std.testing.allocator;

    const databaseType = Database(alphabet.DNA, 8);

    var sequence = try Sequence(alphabet.DNA).init(allocator, "one", "ATCGGTTG");
    defer sequence.deinit();

    var sequences = [_]Sequence(alphabet.DNA){sequence};

    var database = try databaseType.init(allocator, &sequences);
    defer database.deinit();

    var search = try Search(databaseType).init(allocator, database, .{});
    defer search.deinit();
}

//     var extend_align = ea.ExtendAlign(alphabetChosen).init(allocator, .{});
//     defer extend_align.deinit();

//     var banded_align = ba.BandedAlign(alphabetChosen).init(allocator, .{});
//     defer banded_align.deinit();

//     const search_params = SearchParams{};

//     const default_min_hsp_length = 16;
//     const max_hsp_join_distance = 16;

//     var query = try Sequence(alphabetChosen).init(allocator, "test", "AAAGCGCGAGTGAACGCAA");
//     defer query.deinit();

//     const min_hsp_length = std.math.min(default_min_hsp_length, query.data.len / 2);

//     _ = min_hsp_length;
//     _ = max_hsp_join_distance;
//     _ = search_params;

//     var highscores = try Highscores.init(allocator, search_params.max_accepts + search_params.max_rejects);
//     defer highscores.deinit();

//     var kmer_it = bio.kmer.Iterator(kmerInfo).init(query.data);

//     var kmers = try allocator.alloc(kmerInfo.KmerType, kmer_it.num_total());
//     defer allocator.free(kmers);

//     var hits_by_seq = try allocator.alloc(usize, db.sequences.len);
//     defer allocator.free(hits_by_seq);

//     // reset hit counter to 0s
//     std.mem.set(usize, hits_by_seq, 0);

//     var unique_check = try std.DynamicBitSet.initEmpty(allocator, kmerInfo.MaxKmers);
//     defer unique_check.deinit();
