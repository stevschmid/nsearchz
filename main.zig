const std = @import("std");
const print = std.debug.print;
const ArrayList = std.ArrayList;
const Allocator = std.mem.Allocator;

const utils = @import("utils.zig");
const ea = @import("extend_align.zig");
const ba = @import("banded_align.zig");

const Cigar = @import("cigar.zig").Cigar;
const CigarOp = @import("cigar.zig").CigarOp;

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
    var extend_align = ea.ExtendAlign(alphabetChosen).init(allocator, .{});
    defer extend_align.deinit();

    var banded_align = ba.BandedAlign(alphabetChosen).init(allocator, .{});
    defer banded_align.deinit();

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

        // Determine SPS
        sps.clearRetainingCapacity();

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
                        kmers[cursor] == other_kmers[other_cursor])
                    {
                        cursor += 1;
                        other_cursor += 1;
                        length += 1;
                    }

                    // add sps
                    try sps.append(.{ .start_one = pos, .end_one = cursor - 1, .start_two = other_pos, .end_two = other_cursor - 1 });
                }
            }
        }

        // Use the SPS to find all HSP
        // Sort by length
        // Try to find best chain
        // Fill space between with banded align
        var cigar_store = std.AutoHashMap(HSP, Cigar).init(allocator);
        defer {
            var it = cigar_store.iterator();
            while (it.next()) |kv| {
                kv.value_ptr.*.deinit();
            }
            cigar_store.deinit();
        }

        var hsps = std.PriorityDequeue(HSP, void, HSP.lessThanScore).init(allocator, {});
        defer hsps.deinit();

        for (sps.items) |sp| {
            var left_cigar = Cigar.init(allocator);
            defer left_cigar.deinit();

            var right_cigar = Cigar.init(allocator);
            defer right_cigar.deinit();

            var left_result = try extend_align.process(query, seq, .backward, sp.start_one, sp.start_two, &left_cigar);
            var right_result = try extend_align.process(query, seq, .forward, sp.end_one + 1, sp.end_two + 1, &right_cigar);

            var hsp = HSP{ .start_one = left_result.pos_one, .start_two = left_result.pos_two, .end_one = right_result.pos_one, .end_two = right_result.pos_two };

            if (hsp.length() < min_hsp_length)
                continue;

            var pos_one = hsp.start_one;
            var pos_two = hsp.start_two;

            // Construct full hsp (spaced seeds so we cannot assume full match)
            var score: i32 = undefined;

            score = 0;

            var cigar = Cigar.init(allocator);

            score += left_result.score;
            try cigar.appendOther(left_cigar);

            while (pos_one <= hsp.end_one and pos_two <= hsp.end_two) {
                const letter_one = query.data[pos_one];
                const letter_two = seq.data[pos_two];
                const op: CigarOp = if (alphabetChosen.match(letter_one, letter_two)) .match else .mismatch;

                score += alphabetChosen.score(letter_one, letter_two);
                try cigar.add(op);

                pos_one += 1;
                pos_two += 1;
            }

            score += right_result.score;
            try cigar.appendOther(right_cigar);

            // save score
            hsp.score = score;

            var cigar_str = try cigar.toStringAlloc(allocator);
            defer allocator.free(cigar_str);
            std.debug.print("Save {}: {s}\n", .{ hsp.score, cigar_str });

            // save cigar (store will free)
            try cigar_store.put(hsp, cigar);

            // save final hsp
            try hsps.add(hsp);
        }

        // Greedy join HSPs if close
        var chain = std.PriorityQueue(HSP, void, HSP.lessThanPos).init(allocator, {});
        defer chain.deinit();

        // Go through HSP (highest first)
        while (hsps.removeMaxOrNull()) |hsp| {
            var it = chain.iterator();

            var is_overlapping: bool = while (it.next()) |existing_hsp| {
                if (hsp.is_overlapping(existing_hsp))
                    break true;
            } else false;

            if (is_overlapping)
                continue;

            // check if hsp joinable
            it.reset();
            var is_joinable: bool = while (it.next()) |existing_hsp| {
                if (hsp.distance_to(existing_hsp) <= max_hsp_join_distance)
                    break true;
            } else false;

            if (!is_joinable and chain.len > 0)
                continue;

            try chain.add(hsp);
        }

        std.debug.print("Chain size {}\n", .{chain.len});

        // Banded align between chain
        var accept = false;

        if (chain.len > 0) {
            const first_hsp = &chain.items[0];
            const last_hsp = &chain.items[chain.len - 1];

            var final_cigar = Cigar.init(allocator);
            defer final_cigar.deinit();
            var cigar = Cigar.init(allocator);
            defer cigar.deinit();

            // Align first HSP's start to whole sequences begin
            _ = try banded_align.process(query, seq, .backward, first_hsp.start_one, first_hsp.start_two, null, null, &cigar);
            try final_cigar.appendOther(cigar);

            var cursor: usize = 0;
            while (cursor < chain.len) : (cursor += 1) {
                const current_hsp = &chain.items[cursor];
                const current_cigar = cigar_store.get(current_hsp.*).?;
                try final_cigar.appendOther(current_cigar);

                if (cursor + 1 < chain.len) {
                    const next_hsp = &chain.items[cursor + 1];
                    _ = try banded_align.process(query, seq, .forward, current_hsp.end_one + 1, current_hsp.end_two + 1, next_hsp.start_one, next_hsp.start_two, &cigar);
                    try final_cigar.appendOther(cigar);
                }
            }

            // Align last HSP's end to whole sequences end
            _ = try banded_align.process(query, seq, .forward, last_hsp.end_one + 1, last_hsp.end_two + 1, null, null, &cigar);
            try final_cigar.appendOther(cigar);

            var cigar_str = try final_cigar.toStringAlloc(allocator);
            defer allocator.free(cigar_str);

            std.debug.print("Cigar {s}\n", .{cigar_str});
            std.debug.print("Identity {d:.2}\n\n", .{final_cigar.identity()});

            if (final_cigar.identity() >= search_params.min_identity) {
                accept = true;
            }
        }

        if (accept) {
            std.debug.print("HIT!\n", .{});
        } else {
            std.debug.print("MISS!\n", .{});
        }
    } // each candidate
}
