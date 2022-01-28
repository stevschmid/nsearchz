const std = @import("std");
const print = std.debug.print;
const ArrayList = std.ArrayList;
const Allocator = std.mem.Allocator;

const utils = @import("utils.zig");

const Sequence = @import("sequence.zig").Sequence;
const FastaReader = @import("fasta_reader.zig").FastaReader;
const Database = @import("database.zig").Database;
const Highscores = @import("highscores.zig").Highscores;
const HSP = @import("hsp.zig").HSP;

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
    min_identity: f32 = 0.8,
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
    const search_params = SearchParams{};

    const default_min_hsp_length = 16;
    const max_hsp_join_distance = 16;

    var query = try Sequence(alphabetChosen).init(allocator, "test", "AAAGCGCGAGTGAACGCAA");
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

    var sps = std.ArrayList(HSP).init(allocator);
    defer sps.deinit();

    // For each candidate:
    // - Get HSPs,
    // - Check for good HSP (>= similarity threshold)
    // - Join HSP together
    // - Align
    // - Check similarity

    const top_to_bottom = highscores.top_to_bottom();
    for (top_to_bottom) |candidate| {
        const seq_id = candidate.id;
        const seq = db.sequences[seq_id];

        print("Highscore for sequence {s}: {}\n", .{ seq.identifier, candidate.score });

        const offset = db.kmer_offset_by_seq[seq_id];
        const count = db.kmer_count_by_seq[seq_id];
        const other_kmers = db.kmers[offset..(offset + count)];

        for (kmers) |kmer, pos| {
            for (other_kmers) |other_kmer, other_pos| {
                if (kmer != other_kmer)
                    continue;

                if ((pos == 0) or (other_pos == 0) or kmers[pos - 1] == kmerInfo.AmbiguousKmer or other_kmers[other_pos - 1] == kmerInfo.AmbiguousKmer or (kmers[pos - 1] != other_kmers[other_pos - 1])) {
                    var cursor: usize = pos + 1;
                    var other_cursor: usize = other_pos + 1;
                    var length: usize = kmerInfo.NumLettersPerKmer;

                    while (cursor < kmers.len and
                        other_cursor < other_kmers.len and
                        kmers[cursor] != kmerInfo.AmbiguousKmer and
                        other_kmers[other_cursor] != kmerInfo.AmbiguousKmer and
                        kmers[cursor] == other_kmers[other_cursor]) : ({
                        cursor += 1;
                        other_cursor += 1;
                        length += 1;
                    }) {}

                    // add sps
                    try sps.append(HSP{ .a_start = pos, .a_end = cursor - 1, .b_start = other_pos, .b_end = other_cursor + 1 });
                }
            }
        }

        // Use the SPS to find all HSP
        // Sort by length
        // Try to find best chain
        // Fill space between with banded align

    }
}
