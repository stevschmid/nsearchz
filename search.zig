const std = @import("std");
const alphabet = @import("bio/bio.zig").alphabet;

const Database = @import("database.zig").Database;
const Sequence = @import("sequence.zig").Sequence;
const Highscores = @import("highscores.zig").Highscores;

const HSP = @import("hsp.zig").HSP;
const Cigar = @import("cigar.zig").Cigar;
const CigarOp = @import("cigar.zig").CigarOp;

const bio = @import("bio/bio.zig");
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
        kmers: std.ArrayList(kmerInfo.KmerType),
        hits_by_seq: []usize,

        const AlignPart = struct {
            hsp: HSP,
            score: i32,
            cigar: Cigar,

            fn cmpScoreDesc(context: void, left: AlignPart, right: AlignPart) bool {
                _ = context;
                return (left.score > right.score);
            }

            fn cmpPosAsc(context: void, left: AlignPart, right: AlignPart) bool {
                _ = context;
                return (left.hsp.start_one < right.hsp.start_one and left.hsp.start_two < right.hsp.start_two);
            }
        };

        pub fn init(allocator: std.mem.Allocator, database: DatabaseType, options: SearchOptions) !Self {
            return Self{
                .allocator = allocator,
                .database = database,
                .options = options,
                .extend_align = ea.ExtendAlign(A).init(allocator, .{}),
                .banded_align = ba.BandedAlign(A).init(allocator, .{}),
                .kmers = std.ArrayList(kmerInfo.KmerType).init(allocator),
                .hits_by_seq = try allocator.alloc(usize, database.sequences.len),
            };
        }

        pub fn process(self: *Self, query: Sequence(A)) !void {
            const min_hsp_length = std.math.min(DefaultMinHspLength, query.data.len / 2);
            const max_hsp_join_distance = DefaultMaxHSPJoinDistance;
            _ = min_hsp_length;
            _ = max_hsp_join_distance;

            var highscores = try Highscores.init(self.allocator, self.options.max_accepts + self.options.max_rejects);
            defer highscores.deinit();

            // unique tracking
            var unique_check = try std.DynamicBitSet.initEmpty(self.allocator, kmerInfo.MaxKmers);
            defer unique_check.deinit();

            // prepare kmers for query
            var kmer_it = bio.kmer.Iterator(kmerInfo).init(query.data);
            try self.kmers.resize(kmer_it.num_total());
            const kmers = self.kmers.items;

            // reset hit counter to 0s
            const hits_by_seq = self.hits_by_seq;
            std.mem.set(usize, hits_by_seq, 0);

            // lets check these kmers
            var index: usize = 0;
            while (kmer_it.next()) |kmer| : (index += 1) {
                kmers[index] = kmer;

                if (kmer == kmerInfo.AmbiguousKmer or unique_check.isSet(kmer))
                    continue;

                unique_check.set(kmer);

                const offset = self.database.seq_offset_by_kmer[kmer];
                const count = self.database.seq_count_by_kmer[kmer];

                const seq_indices = self.database.seq_indices[offset..(offset + count)];
                for (seq_indices) |seq_index| {
                    hits_by_seq[seq_index] += 1;
                    highscores.add(seq_index, hits_by_seq[seq_index]);
                }
            }

            var sps = std.ArrayList(HSP).init(self.allocator);
            defer sps.deinit();

            var num_hits: usize = 0;
            var num_rejects: usize = 0;

            const top_to_bottom = highscores.top_to_bottom();
            for (top_to_bottom) |candidate| {
                const seq_id = candidate.id;
                const seq = self.database.sequences[seq_id];

                const offset = self.database.kmer_offset_by_seq[seq_id];
                const count = self.database.kmer_count_by_seq[seq_id];
                const other_kmers = self.database.kmers[offset..(offset + count)];

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

                var align_parts = std.ArrayList(AlignPart).init(self.allocator);
                defer {
                    for (align_parts.items) |*align_part| {
                        align_part.cigar.deinit();
                    }
                    align_parts.deinit();
                }

                for (sps.items) |sp| {
                    var left_cigar = Cigar.init(self.allocator);
                    defer left_cigar.deinit();

                    var right_cigar = Cigar.init(self.allocator);
                    defer right_cigar.deinit();

                    var left_result = try self.extend_align.process(query, seq, .backward, sp.start_one, sp.start_two, &left_cigar);
                    var right_result = try self.extend_align.process(query, seq, .forward, sp.end_one + 1, sp.end_two + 1, &right_cigar);

                    var hsp = HSP{ .start_one = left_result.pos_one, .start_two = left_result.pos_two, .end_one = right_result.pos_one, .end_two = right_result.pos_two };

                    if (hsp.length() < min_hsp_length)
                        continue;

                    var pos_one = hsp.start_one;
                    var pos_two = hsp.start_two;

                    // Construct full hsp (spaced seeds so we cannot assume full match)
                    var score: i32 = 0;
                    score = left_result.score;

                    var cigar = Cigar.init(self.allocator);
                    try cigar.appendOther(left_cigar);

                    // go until we hit start of SP (not HSP)
                    while (pos_one <= sp.end_one and pos_two <= sp.end_two) {
                        const letter_one = query.data[pos_one];
                        const letter_two = seq.data[pos_two];
                        const op: CigarOp = if (A.match(letter_one, letter_two)) .match else .mismatch;

                        score += A.score(letter_one, letter_two);
                        try cigar.add(op);

                        pos_one += 1;
                        pos_two += 1;
                    }

                    score += right_result.score;
                    try cigar.appendOther(right_cigar);

                    const align_part = AlignPart{ .hsp = hsp, .cigar = cigar, .score = score };
                    try align_parts.append(align_part);
                }

                // Sort by score, highest first
                std.sort.sort(AlignPart, align_parts.items, {}, AlignPart.cmpScoreDesc);

                // Greedy join HSPs if close
                var chain = std.ArrayList(AlignPart).init(self.allocator);
                defer chain.deinit();

                // Go through HSP (highest first)
                for (align_parts.items) |part| {
                    // check if overlapping
                    var is_overlapping: bool = for (chain.items) |chain_part| {
                        if (part.hsp.is_overlapping(chain_part.hsp))
                            break true;
                    } else false;

                    if (is_overlapping)
                        continue;

                    // check if hsp joinable
                    var is_joinable: bool = for (chain.items) |chain_part| {
                        if (part.hsp.distance_to(chain_part.hsp) <= max_hsp_join_distance)
                            break true;
                    } else false;

                    if (!is_joinable and chain.items.len > 0)
                        continue;

                    try chain.append(part);
                }

                // Sort by Pos
                std.sort.sort(AlignPart, chain.items, {}, AlignPart.cmpPosAsc);

                // Banded align between chain
                var accept = false;

                if (chain.items.len > 0) {
                    const first_part = &chain.items[0];
                    const last_part = &chain.items[chain.items.len - 1];

                    var cigar = Cigar.init(self.allocator);
                    defer cigar.deinit();

                    var final_cigar = Cigar.init(self.allocator);
                    defer final_cigar.deinit();

                    // Align first HSP's start to whole sequences begin
                    _ = try self.banded_align.process(query, seq, .backward, first_part.hsp.start_one, first_part.hsp.start_two, null, null, &cigar);
                    try final_cigar.appendOther(cigar);

                    for (chain.items) |part, part_index| {
                        try final_cigar.appendOther(part.cigar);

                        if (part_index + 1 < chain.items.len) {
                            const next_part = &chain.items[part_index + 1];
                            _ = try self.banded_align.process(query, seq, .forward, part.hsp.end_one + 1, part.hsp.end_two + 1, next_part.hsp.start_one, next_part.hsp.start_two, &cigar);
                            try final_cigar.appendOther(cigar);
                        }
                    }

                    // Align last HSP's end to whole sequences end
                    _ = try self.banded_align.process(query, seq, .forward, last_part.hsp.end_one + 1, last_part.hsp.end_two + 1, null, null, &cigar);
                    try final_cigar.appendOther(cigar);

                    var cigar_str = try final_cigar.toStringAlloc(self.allocator);
                    defer self.allocator.free(cigar_str);

                    if (final_cigar.identity() >= self.options.min_identity) {
                        accept = true;
                    }
                }

                if (accept) {
                    num_hits += 1;
                } else {
                    num_rejects += 1;
                }

                if (num_hits >= self.options.max_accepts or num_rejects >= self.options.max_rejects)
                    break;
            } // each candidate
        }

        pub fn deinit(self: *Self) void {
            self.extend_align.deinit();
            self.banded_align.deinit();
            self.kmers.deinit();
            self.allocator.free(self.hits_by_seq);
        }
    };
}

test "check" {
    const allocator = std.testing.allocator;

    const sequences = [_]Sequence(alphabet.DNA){
        .{ .identifier = "DB1", .data = "ATCGTGAGACGATGCAAAAAATTGAGA" },
        .{ .identifier = "DB2", .data = "GTCCGACGCAATAAACTATATGGGG" },
    };
    const query = Sequence(alphabet.DNA){ .identifier = "Query", .data = "GGTGAGACGACGCAATAAATTGAGA" };

    const databaseType = Database(alphabet.DNA, 8);
    var database = try databaseType.init(allocator, &sequences);
    defer database.deinit();

    var search = try Search(databaseType).init(allocator, database, .{});
    try search.process(query);
    defer search.deinit();
}
