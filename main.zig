const std = @import("std");
const print = std.debug.print;
const ArrayList = std.ArrayList;
const Allocator = std.mem.Allocator;

const utils = @import("utils.zig");

const Sequence = @import("sequence.zig").Sequence;
const FastaReader = @import("fasta_reader.zig").FastaReader;
const Database = @import("database.zig").Database;
const Highscores = @import("highscores.zig").Highscores;

const bio = @import("bio/bio.zig");
const alphabet = bio.alphabet;

const alphabetChosen = alphabet.DNA;
const kmerInfo = bio.kmer.KmerInfo(alphabetChosen, 8);

const Strand = enum {
    Plus,
    Minus,
    Both,
};

const SearchParams = struct {
    max_accepts: u32 = 1,
    max_rejects: u32 = 16,
    min_identity: f32 = 0.75,
    strand: Strand = Strand.Plus,
};

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    var bench_start = std.time.milliTimestamp();
    var reader = FastaReader(alphabetChosen).init(allocator);
    defer reader.deinit();

    var arg_it = std.process.args();

    // skip my own exe name
    _ = arg_it.skip();

    const file = (try arg_it.next(allocator) orelse {
        print("Expected first argument to be path to input file\n", .{});
        return error.InvalidArgs;
    });
    defer allocator.free(file);
    try reader.readFile(file);

    print("Reading took {}ms\n", .{std.time.milliTimestamp() - bench_start});

    bench_start = std.time.milliTimestamp();
    var db = try Database(alphabet.DNA, 8).init(allocator, reader.sequences.items);
    defer db.deinit();

    print("Indexing took {}ms\n", .{std.time.milliTimestamp() - bench_start});

    // search 
    const search_params = SearchParams {};

    const default_min_hsp_length = 16;
    const max_hsp_join_distance = 16;

    var query = try Sequence(alphabetChosen).init(allocator, "test", "CGATTACCGTTGATTTCA");
    defer query.deinit();

    const min_hsp_length = std.math.min(default_min_hsp_length, query.data.len / 2);

    _ = min_hsp_length;
    _ = max_hsp_join_distance;
    _ = search_params;

    var highscores = try Highscores.init(allocator, search_params.max_accepts + search_params.max_rejects);
    defer highscores.deinit();

    var kmer_it = bio.kmer.Iterator(kmerInfo).init(query.data);

    var kmers = try allocator.alloc(kmerInfo.KmerType, kmer_it.num_total());
    defer allocator.free(kmers);

    var hits_by_seq = try allocator.alloc(usize, db.sequences.len);
    defer allocator.free(hits_by_seq);

    // reset hit counter to 0s
    std.mem.set(usize, hits_by_seq, 0);

    var unique_check = try std.DynamicBitSet.initEmpty(allocator, kmerInfo.MaxKmers);
    defer unique_check.deinit();

    var index: usize = 0;
    while (kmer_it.next()) |kmer| : (index += 1) {
        kmers[index] = kmer;

        if (kmer == kmerInfo.AmbiguousKmer or unique_check.isSet(kmer))
            continue;

        unique_check.set(kmer);

        const offset = db.seq_offset_by_kmer[kmer];
        const count = db.seq_count_by_kmer[kmer];

        const seq_indices = db.seq_indices[offset..(offset + count)];
        for (seq_indices) |seq_index| {
            // highscores.add(seq_index, 
            hits_by_seq[seq_index] += 1;
            highscores.add(seq_index, hits_by_seq[seq_index]);
        }
    }

    for (highscores.result()) |highscore_entry| {
        print("Highscore {} -> {}\n", .{highscore_entry.id, highscore_entry.score});
    }
}
