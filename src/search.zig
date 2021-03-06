const std = @import("std");
const alphabet = @import("bio/bio.zig").alphabet;

const Database = @import("database.zig").Database;
const Sequence = @import("sequence.zig").Sequence;
const SequenceList = @import("sequence.zig").SequenceList;
const Highscores = @import("highscores.zig").Highscores;

const HSP = @import("hsp.zig").HSP;
const Cigar = @import("cigar.zig").Cigar;
const CigarOp = @import("cigar.zig").CigarOp;

const bio = @import("bio/bio.zig");
const ea = @import("extend_align.zig");
const ba = @import("banded_align.zig");

const utils = @import("utils.zig");

pub const Strand = enum {
    plus,
    minus,
    both,
};

pub const SearchOptions = struct {
    max_accepts: u32 = 1,
    max_rejects: u32 = 16,
    min_identity: f32 = 0.75,
    strand: Strand = .both,
};

pub fn SearchHit(comptime A: type) type {
    return struct {
        const Self = @This();

        allocator: std.mem.Allocator,
        target: Sequence(A),
        cigar: Cigar,
        reverse_match: bool = false,

        pub fn init(allocator: std.mem.Allocator, target: Sequence(A), cigar: Cigar, reverse_match: bool) !Self {
            return Self{
                .allocator = allocator,
                .target = try target.clone(),
                .cigar = try cigar.clone(),
                .reverse_match = reverse_match,
            };
        }

        pub fn clone(self: Self) !Self {
            return try Self.init(self.allocator, self.target, self.cigar, self.reverse_match);
        }

        pub fn deinit(self: *Self) void {
            self.target.deinit();
            self.cigar.deinit();
        }
    };
}

pub fn SearchHitList(comptime A: type) type {
    return utils.ArrayListDeinitWrapper(SearchHit(A));
}

const ArrayPartList = utils.ArrayListDeinitWrapper(AlignPart);

const AlignPart = struct {
    const Self = @This();

    allocator: std.mem.Allocator,
    hsp: HSP,
    score: i32,
    cigar: Cigar,

    pub fn init(allocator: std.mem.Allocator) Self {
        return Self{
            .allocator = allocator,
            .cigar = Cigar.init(allocator),
            .hsp = undefined,
            .score = undefined,
        };
    }

    pub fn deinit(self: *Self) void {
        self.cigar.deinit();
    }

    fn cmpScoreDesc(context: void, left: AlignPart, right: AlignPart) bool {
        _ = context;
        return (left.score > right.score);
    }

    fn cmpPosAsc(context: void, left: AlignPart, right: AlignPart) bool {
        _ = context;
        return (left.hsp.start_one < right.hsp.start_one and left.hsp.start_two < right.hsp.start_two);
    }
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
        database: *DatabaseType,
        extend_align: ea.ExtendAlign(A),
        banded_align: ba.BandedAlign(A),
        kmers: std.ArrayList(kmerInfo.KmerType),
        hits_by_seq: []usize,

        pub fn init(allocator: std.mem.Allocator, database: *DatabaseType, options: SearchOptions) !Self {
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

        pub fn search(self: *Self, query: Sequence(A), hits: *SearchHitList(A)) !void {
            hits.list.clearRetainingCapacity();

            if (A.SupportsStrands) {
                if (self.options.strand == .plus or self.options.strand == .both) {
                    try self.do(query, hits, false);
                }

                if (self.options.strand == .minus or self.options.strand == .both) {
                    var reverse_complement = try query.clone();
                    defer reverse_complement.deinit();

                    reverse_complement.reverse();
                    reverse_complement.complement();

                    try self.do(reverse_complement, hits, true);
                }
            } else {
                try self.do(query, hits, false);
            }
        }

        fn do(self: *Self, query: Sequence(A), hits: *SearchHitList(A), reverse_match: bool) !void {
            const min_hsp_length = std.math.min(DefaultMinHspLength, query.data.len / 2);
            const max_hsp_join_distance = DefaultMaxHSPJoinDistance;

            // highscores
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

            var left_cigar = Cigar.init(self.allocator);
            defer left_cigar.deinit();

            var right_cigar = Cigar.init(self.allocator);
            defer right_cigar.deinit();

            var cigar = Cigar.init(self.allocator);
            defer cigar.deinit();

            var final_cigar = Cigar.init(self.allocator);
            defer final_cigar.deinit();

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

                var align_parts = ArrayPartList.init(self.allocator);
                defer align_parts.deinit();

                for (sps.items) |sp| {
                    // TODO: can we use already aligned parts to check out if this seed is part of the an already extended seed, so we can save
                    var left_result = try self.extend_align.process(query, seq, .backward, sp.start_one, sp.start_two, &left_cigar);
                    var right_result = try self.extend_align.process(query, seq, .forward, sp.end_one + 1, sp.end_two + 1, &right_cigar);

                    var align_part = AlignPart.init(self.allocator);
                    align_part.hsp = .{ .start_one = left_result.pos_one, .start_two = left_result.pos_two, .end_one = right_result.pos_one, .end_two = right_result.pos_two };

                    if (align_part.hsp.length() < min_hsp_length)
                        continue;

                    var pos_one = sp.start_one;
                    var pos_two = sp.start_two;

                    try align_part.cigar.appendOther(left_cigar);

                    // Construct full hsp (spaced seeds so we cannot assume full match)
                    align_part.score = left_result.score;

                    // go until we hit start of SP (not HSP)
                    while (pos_one <= sp.end_one and pos_two <= sp.end_two) {
                        const letter_one = query.data[pos_one];
                        const letter_two = seq.data[pos_two];
                        const op: CigarOp = if (A.match(letter_one, letter_two)) .match else .mismatch;

                        align_part.score += A.score(letter_one, letter_two);
                        try align_part.cigar.add(op);

                        pos_one += 1;
                        pos_two += 1;
                    }

                    align_part.score += right_result.score;
                    try align_part.cigar.appendOther(right_cigar);

                    try align_parts.list.append(align_part);
                }

                // Sort by score, highest first
                std.sort.sort(AlignPart, align_parts.list.items, {}, AlignPart.cmpScoreDesc);

                // Greedy join HSPs if close
                var chain = std.ArrayList(AlignPart).init(self.allocator);
                defer chain.deinit();

                // Go through HSP (highest first)
                for (align_parts.list.items) |part| {
                    // check if overlapping
                    var is_overlapping: bool = for (chain.items) |chain_part| {
                        if (part.hsp.is_overlapping(chain_part.hsp))
                            break true;
                    } else false;

                    if (is_overlapping)
                        continue;

                    // check if hsp joinable
                    var is_joinable: bool = for (chain.items) |chain_part| {
                        const distance_satisified = part.hsp.distance_to(chain_part.hsp) <= max_hsp_join_distance;
                        const downstream_satisfied = part.hsp.is_downstream_of(chain_part.hsp);

                        if (distance_satisified and downstream_satisfied)
                            break true;
                    } else false;

                    if (!is_joinable and chain.items.len > 0)
                        continue;

                    try chain.append(part);
                }

                // Sort by Pos
                std.sort.sort(AlignPart, chain.items, {}, AlignPart.cmpPosAsc);

                // Banded align between chain
                if (chain.items.len > 0) {
                    const first_part = &chain.items[0];
                    const last_part = &chain.items[chain.items.len - 1];

                    final_cigar.clear();

                    // Align first HSP's start to whole sequences begin
                    _ = try self.banded_align.process(query, seq, .backward, first_part.hsp.start_one, first_part.hsp.start_two, null, null, &cigar);
                    try final_cigar.appendOther(cigar);

                    index = 0;
                    while (index + 1 < chain.items.len) : (index += 1) {
                        const part = chain.items[index];
                        const next_part = chain.items[index + 1];

                        try final_cigar.appendOther(part.cigar);
                        _ = try self.banded_align.process(
                            query,
                            seq,
                            .forward,
                            part.hsp.end_one + 1,
                            part.hsp.end_two + 1,
                            next_part.hsp.start_one,
                            next_part.hsp.start_two,
                            &cigar,
                        );
                        try final_cigar.appendOther(cigar);
                    }

                    // Align last HSP's end to whole sequences end
                    try final_cigar.appendOther(last_part.cigar);
                    _ = try self.banded_align.process(query, seq, .forward, last_part.hsp.end_one + 1, last_part.hsp.end_two + 1, null, null, &cigar);
                    try final_cigar.appendOther(cigar);

                    var accept = (final_cigar.identity() >= self.options.min_identity);
                    if (accept) {
                        try hits.list.append(try SearchHit(A).init(hits.allocator, seq, final_cigar, reverse_match));

                        num_hits += 1;
                    } else {
                        num_rejects += 1;
                    }
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

const databaseType = Database(alphabet.DNA, 8);

test "check" {
    const allocator = std.testing.allocator;

    var sequences = SequenceList(alphabet.DNA).init(allocator);
    defer sequences.deinit();
    try sequences.list.append(try Sequence(alphabet.DNA).init(allocator, "DB1", "ATCGTGAGACGATGCAAAAAATTGAGA"));
    try sequences.list.append(try Sequence(alphabet.DNA).init(allocator, "DB2", "GTCCGACGCAATAAACTATATGGGG"));

    var query = try Sequence(alphabet.DNA).init(allocator, "Query", "GGTGAGACGACGCAATAAATTGAGA");
    defer query.deinit();

    var database = try databaseType.init(allocator, sequences.list.toOwnedSlice(), null);
    defer database.deinit();

    {
        var search = try Search(databaseType).init(allocator, &database, .{});
        defer search.deinit();

        var hits = SearchHitList(alphabet.DNA).init(allocator);
        defer hits.deinit();

        try search.search(query, &hits);

        try std.testing.expectEqual(@as(usize, 1), hits.list.items.len);
        try std.testing.expectEqualStrings("DB1", hits.list.items[0].target.identifier);
    }

    // try accepts 2, other sequence is still low
    {
        var search = try Search(databaseType).init(allocator, &database, .{ .max_accepts = 2 });
        defer search.deinit();

        var hits = SearchHitList(alphabet.DNA).init(allocator);
        defer hits.deinit();

        try search.search(query, &hits);

        // still 1
        try std.testing.expectEqual(@as(usize, 1), hits.list.items.len);
    }

    // accept two, but lower identity threshold
    {
        var search = try Search(databaseType).init(allocator, &database, .{ .max_accepts = 2, .min_identity = 0.5 });
        defer search.deinit();

        var hits = SearchHitList(alphabet.DNA).init(allocator);
        defer hits.deinit();

        try search.search(query, &hits);

        // now 2
        try std.testing.expectEqual(@as(usize, 2), hits.list.items.len);
    }
}

test "search multiple" {
    const allocator = std.testing.allocator;

    var sequences = SequenceList(alphabet.DNA).init(allocator);
    defer sequences.deinit();
    try sequences.list.append(try Sequence(alphabet.DNA).init(allocator, "DB1", "ATCGTGAGACGATGCAAAAAATTGAGA"));
    try sequences.list.append(try Sequence(alphabet.DNA).init(allocator, "DB2", "GTCCGACGCAATAAACTATATGGGG"));

    var database = try databaseType.init(allocator, sequences.list.toOwnedSlice(), null);
    defer database.deinit();

    var search = try Search(databaseType).init(allocator, &database, .{});
    defer search.deinit();

    var query1 = try Sequence(alphabet.DNA).init(allocator, "Query 1 ", "TGAGACGATGCAAA");
    defer query1.deinit();

    var hits1 = SearchHitList(alphabet.DNA).init(allocator);
    defer hits1.deinit();

    try search.search(query1, &hits1);
    try std.testing.expectEqual(@as(usize, 1), hits1.list.items.len);

    var query2 = try Sequence(alphabet.DNA).init(allocator, "Query 1 ", "CGTTATATTCGGAGACCTAT");
    defer query2.deinit();

    var hits2 = SearchHitList(alphabet.DNA).init(allocator);
    defer hits2.deinit();

    try search.search(query2, &hits2);
    try std.testing.expectEqual(@as(usize, 0), hits2.list.items.len);
}

test "ultrasequence" {
    const allocator = std.testing.allocator;

    var sequences = SequenceList(alphabet.DNA).init(allocator);
    defer sequences.deinit();

    try sequences.list.append(try Sequence(alphabet.DNA).init(allocator, "RF01960;SSU_rRNA_eukarya;AADB02014478.1/16479-14622   9606:Homo sapiens (human)", "UACAUGGUUGAUCCUACCAGAACGAUAUGCUUGUCUCAAAGAUUAAGCCAUACAUGUCUAAGUACGCAGGGCCGGUACAGUGAAACUGCGAAUGGCUCAUUAAAUCAGUUAUGGUUCUUUUGUUUGCUCGCUCCUCUCCUACUUGGAAAACUGUGGUAAUUCUAGAGCUAAUACAUGCCAAAGGGCGCUGACCCCCUUCGCGGGGAAGAUGCGUGCAUUUAGCAGAUCAAAACCAACCCAGUCAGCCCCUCUCCAGCCCCGGCUGGAGGGUCAGGUGCCACUGGCUUUGGUGACUCUAGAUAACCUCAGGCCAAUUGCACGCCCCCAGUGGCAGCGAUGACCCAUUGUAAAGUCUGCCCUAUCAACUUUCGAUGGUAGUCGCUGUGCCUACCAUGGUGACCACGGGUGACAAGGAAUCAGGGUUCGAUUCUGGAGAGGGAGCCUGAGAAAUGGCUACCACAUUCAAGGAAGGCAGCAGGCAUGCAAAUUACCCACUCCCUACUCGGGGAGGUAGUGAUGAAAAAUAACAAUACAGGACUCUUUCGAGGCCCUGUAAUUGGAAUGAGUCCACUUUAAAUCCUUUAACGAGGAUCCAUUGGAGGGCAAGUCUGGUGCCAGCAGCUGCAGUAAUUCCAGCUCCAAUAGCGUAUAUUAAAGUUGCUGUAGUUAAAAAGCUCGUAGUUGGAUCUUGGGAGCGGGCGGGCGGGCGAGCCAUGGCCCGUCCCCGCCCCUUGCCUCUUGGCGCCCCCUCGAUGUUCUUAGCUGAGUAUCCCGUGGCCCUAAGCUUUUACUUUGAAAAAGUUAGAGUUUUCAAAGCAGGCCCGGGCCGCCUAGAUACCGCAGCUAGGAAUGAUGGAAUAGAACCGCGGUUCUAUUUUGUUGGUUUUUGGAACUGAGGCCAUGAUUAAGAAGGACGGCUGGGGGCAUUCGUAUUGCGCCACUACAGGUGAAAUUCUUGGACCGGCGCAGUUCGGACCAGAGCGAAAGCAUUUGCCAAGAAUGUUUUCAUUAACCAAGAACGAAAGUCGGAGGUUCGAAGACGAUCAGAUACCAUCGUAGUUCCGACCAUAAAAGAUGCCGACUGGCGAUGCGUGCGGCAGUGUUAUUUUCAUGACCCACUGGGCAGCUUCCAAGAAACCAAAGUCUUUGGGUUUUUGGGUUCCCGGGGGGAGUAUGGUUGCAAAGCUGAAACUUAAAGGAAUUGGUGGAAGGGCACCACCAGGAGUGGAGCCUGUGGCUUAAUUUGACUCAACACGGGAAACCUCACCCGGCCCAGACACAGACAGGAUUGACAGACUGAUAGCUCCUUCUCAAUUCCAUGGGUGGUGGUGCAUGGCCUUAGCUGGUGGAGUUAUUUGUCUGGUUAAUUCCGAUAACGAACGAGACUCUGGCAUGCUAACUAAUUACGUGACCCACAAGCGGUCGGUGUCCCCCAACUUCUUAGAUGGAUAAGUGGCAUUCAGCCACCCGAGAUUGAGUAAUAACAGGUCUGUGAUGCCCUUAGAUGUCCGGGGCUGCACGUGCGCUACACUGACUGGCUCAGCGUGUGCCUACCAUACGUGGCAGACGUGGGUAACCUGUUGAAUACCAUUGGUGAUGGUGAUCGGGGAUUGCAAUUAUUCCCCAUGAACGAGGAAUUCCCAGUAAGUGUGGGUCAUAGGCUUGUGUGAAGUCCCUACCCUUUGUACACACCGCCAGUUGCUACUACUGAUUGGAUGGUUUAGUGAGGCCCUCAGAGCGGCCCCGCCGGGUCAGCCCACUGCCGUGGCCGAGCACUGAGAAGACAGUCGAACUUGACUAUCUAGAGGAAGUAAAAGUCGUAACAAAGUUUCUGUAGGUGAACCUGUGGAGGGAUCAU"));

    var query = try Sequence(alphabet.DNA).init(allocator, "RF01960;SSU_rRNA_eukarya;CABZ01060450.1/90976-89252   7955:Danio rerio (zebrafish) ", "CUUUUGGCCAAUUUAACCUUGGCACUUUAAUAUUUCUGAACAAGAUACACUAGUACAUGUGGGCUCAUUAGAGAGCUAAGAAUGUGUUUUUUAAUGAUGUUGGUUAUUUAAAAAUAGCUUAAAUCUUCCCUGAGUAAUCAAAAAAAAAAAAGCUACCAUCAAACAGGCACUCCAAGCACAUGCAUAAUAAAAAAAACAGUCAAAUUUUAAUAGUAUGACCAACAUGUGUGUGUUUGUUUUCAUUUUUUAUACUCUAAAGUGUCAAAGUACUGAGUGCGCUGACUUUUGGCCAAUUUAACCUUGGCACUUUAAUAUUUCUAAACAAGAUACACUAGUACAUGUGGCCUCAUUAGAAAGCUAAGAAUUUAAUGAUGUUGGUUACUUAAAAAUAGCUUAAAUCUUCACUGAUUUAUCAAAAAAAUUCACGCUUUUAAACAGACACUGAAUAAGAACAAUACAGGUCUCUUUCGAGGCCCUGUAACUGGAAUGAGCGUAUCCUAAACCCAUGGGCGAGGACCCAUUAGAGGGCAAGUCUGGUGCCAGCAGCCGCGGUAAUUCCAGCUCCAAUAGCGUAUAUCAAAGUUGCAGCAGUUAGAAAGCUCGUAGUUGGAUUUUGGGAGUGGGAUGGCGGUCCGCUGCCCUUAGCUGGGUGUCCGGUACCUCGGGGCCCAGAGCGUUUACUUUGAAAAAAUUUGAGUGUUCAAAGCAGGCCGGCCAGCCGCCGCUGAAUACCGCAGCUAGGAAUAAUGGAAUAGGACUCCGGUUCUAUUUUCUGGAACCUGGGGCCAUGAUUGAGAGGGACGGCCGGGGGCAUUCGUAUUGCGCCGCUAGAGGUGAAAUUCUUGGACCGGCACAAGACGCACGAAAGCGAAAUCAUUUGCCAAGAAUGUUUUCAUUAAUCAAGAACGAAAGUCGGAGGUUCGAAGACGAUCAGAUACCGUCAUAGUUCCGACCGUAAACAAUGCCGACCCGCGAUCCGGCGGCGUUAUUCCCAUGACCCGCCGGGCAGCGUGCGGGAAACCAUGAGUCUUUGGGUUCUGGGGGGAAUAUGGUUGCAAGGCUGAAACUUGAAGAAAUUGACGGAAGGGCAGCACCAGGAGUGGAGCCUGCGGCUUAAUUUGACUCAACACAGGAAACCUCACCCGGCCUGGACACGGAAAGGAUUGACAGAUUGAUAGCUCUUUCUCGAUUCUGUGGGUGGUGGUGCAUGGCCAUUCUAAGUUAGUGGACUGGUUCAUUCCGAUAACGAACGAGACUCCGGCAUGGUGAGUUAUCAAAAAAGCUACCAUCAAACAGGCACUCCAAGCACAUGUGCAAUAAAAUAAUUGUCAAAUUUGAAUAGUAUGACCACAUGUGUGUAUUUAUUUUCAUUUUUAUACUCUAAAGUGUCAAAGUACUCAGUGAGAUGCAGAAAACCAACAACUUUUAUCCAUUUAUCCAAUAUUUCUACUUCUGAAGAAGACAAAUGUAUGUUAAAUUUGGGUUAUCUGUCAAGUUGAAUUUUGAUUAUUUCAACAGUGAUGGUGUGUUGACGUUUUGACAAUUUAACCUAAAUUUUCCUGAGCUUUCCAUCAAUGUAUGGCAUGAAGCUGCGGGAGAUUUGUUGUUCCUGUCUCUUAAAAUCCCUGUUCAAAUUGAGCCAGUUUUCAGUGAUUUGCAGUGGGCCAGUUAUUGGACGGUAUUAAAAAAUGCACAACUUUGGUUGUGCUUCUGCAAAAAUUAGGU");
    defer query.deinit();

    var database = try databaseType.init(allocator, sequences.list.toOwnedSlice(), null);
    defer database.deinit();

    var search = try Search(databaseType).init(allocator, &database, .{ .min_identity = 0.5 });
    defer search.deinit();

    var hits = SearchHitList(alphabet.DNA).init(allocator);
    defer hits.deinit();

    try search.search(query, &hits);
    try std.testing.expectEqual(@as(usize, 1), hits.list.items.len);

    const hit = hits.list.items[0];

    // expect nsearch output
    try std.testing.expectApproxEqAbs(@as(f32, 0.53), hit.cigar.identity(), 0.01);
}

const databaseContents =
    \\>RF00807;mir-314;AFFE01007792.1/82767-82854   42026:Drosophila bipectinata
    \\UCGUAACUUGUGUGGCUUCGAAUGUACCUAGUUGAGGAAAAAUCAGUUUG
    \\GAUUUUGUUACCUCUGGUAUUCGAGCCAAUAAGUUCGG
    \\>RF00807;mir-314;AAFS01000446.1/64778-64866   46245:Drosophila pseudoobscura pseudoobscura
    \\UCGUAACUUGUGUGGCUUCGAAUGUACCUAGUUGAGGAAAACUCCGAAAU
    \\GGAUUUUGUUACCUCUGGUAUUCGAGCCAAUAAGUUCGG
    \\>RF00807;mir-314;AANI01017486.1/342740-342830   7244:Drosophila virilis
    \\UCGUAACUUGUGUGGCUUGAAUGUACCUGGUUGAGGAACGAAUUCAACGU
    \\UUGGAUUUUGUUGCCUUUGGUAUUCGAGCCAAUAAGUUCGG
    \\>RF00807;mir-314;AAPU01011627.1/156896-156990   7230:Drosophila mojavensis
    \\UCGUAACUUGUGUGGCUUCGAAUGUACCUCGUCGAGCGAAAAGCGAAUUC
    \\AUUGUUGGAUUUUGUUGCUCUUGGUAUUCGAGCCAAUAAGUUCGG
    \\>RF00752;mir-14;AAWU01029067.1/5309-5242   7176:Culex quinquefasciatus (southern house mosquito)
    \\UGUGGGAGCGAGAUUAAGGCUUGCUGGUUUCACGUUCGAGUAAAGUCAGU
    \\CUUUUUCUCUCUCCUAUU
    \\>RF00752;mir-14;AAGE02012112.1/774-706   7159:Aedes aegypti (yellow fever mosquito)
    \\UGUGGGAGCGAGAUUAAGGCUUGCUGGUCAUUUAUUACACUCGAAGUCAG
    \\UCUUUUUCUCUCUCCUAUU
    \\>RF00752;mir-14;AANI01011011.1/11101-11163   7244:Drosophila virilis
    \\UGUGGGAGCGAGACGGGGACUCACUGUGCUUUUUAUAUAGUCAGUCUUUU
    \\UCUCUCUCCUAUA
    \\>RF00715;mir-383;AAMC01036319.1/13960-13888   8364:Xenopus (Silurana) tropicalis (western clawed frog)
    \\CUCCUCAGAUCAGAAGGUGAUUGUGGCUUUUAGUAGAUAUUAAGCAGCCA
    \\CAGCACUGCCUGGUCAGAAAGAG
    \\>RF00715;mir-383;AAQR03137803.1/1757-1829   30611:Otolemur garnettii (small-eared galago)
    \\CUCCUCAGAUCAGAAGGUGAUUGUGGCUUUGGGUGCAUGGUUAUAAGCCA
    \\CAGCACUGCCUGGUCAGAAAGAG
    \\>RF00715;mir-383;AFEY01400405.1/5982-5910   9305:Sarcophilus harrisii (Tasmanian devil)
    \\CUCCUCAGAUCAGAAGGUGAUUGUGGCUUUGGGCAGACAUGGAACAGCCA
    \\CAUCACUGGCUGGUCAGAAAGAG
    \\>RF01157;sn1185;DQ789405.1/1-66   6238:Caenorhabditis briggsae
    \\AUCGGUGAUGUGAUAUCCAGUUCUGCUACUGAAGCGUUGUGAAGAUUAAC
    \\UUUCCCCGUCUGAGAU
    \\>RF01157;sn1185;AEHI01092347.1/1008-1073   860376:Caenorhabditis angaria
    \\ACUGAUGAUGUUAACUCCAGUUCUGCUACUGAAUGAAUGUGACGAUAUUC
    \\UUUCCCCGACUGAGGU
    \\>RF01157;sn1185;ABLE03029241.1/4849-4913   281687:Caenorhabditis japonica
    \\AUUGAUGAUGUUCAUCCAGUUCUGCUACUGAAUCAGUGUGAAGAUAUUCU
    \\UUCCCCGACUGAGAU
    \\>RF01885;HSR-omega_1;AFPP01029324.1/15839-15914   1041015:Drosophila rhopaloa
    \\ACCACCUAACCAAGCAAUAUGUAUUUCUUUCUCUAAACUUUAUAGUUGGG
    \\CGUUGAAAGUUGAUAUCGAUCCGUGA
    \\>RF01885;HSR-omega_1;AFFH01007186.1/256640-256716   30033:Drosophila kikkawai
    \\AUCACUUAACCAGCAAUAUGUAUUUCUUUCUCUAAACUUUAUAGUUGGGC
    \\GUUGAAAGUUGAUACGCGAACGUGAAA
    \\>RF01885;HSR-omega_1;AAPU01011178.1/205842-205767   7230:Drosophila mojavensis
    \\ACACGUUAACCAAGCAUUAUGUAUUUCUUUCUCUAAACUUUAUAGUUGGG
    \\CGUUGAAAGUUGAUACGCGAUCGAAC
    \\
;

test "strand support" {
    const FastaReader = @import("io/fasta_reader.zig").FastaReader;

    const allocator = std.testing.allocator;

    var db_source = std.io.StreamSource{ .const_buffer = std.io.fixedBufferStream(databaseContents) };

    var db_reader = FastaReader(alphabet.DNA).init(allocator, &db_source);
    defer db_reader.deinit();

    var sequences = SequenceList(alphabet.DNA).init(allocator);
    defer sequences.deinit();

    try sequences.list.append(try Sequence(alphabet.DNA).init(allocator, "DB2", "GTCCGACGCAATAAACTATATGGGG"));

    while (try db_reader.next()) |sequence| {
        try sequences.list.append(sequence);
    }

    var database = try databaseType.init(allocator, sequences.list.toOwnedSlice(), null);
    defer database.deinit();

    var query = try Sequence(alphabet.DNA).init(
        allocator,
        "RF00807;mir-314;AAPT01020574.1/773332-773257   7222:Drosophila grimshawi",
        "UCGUAACUUGUGUGGCUUCGAAUGUACCUGGCUAAGGAAAGUUGGAUUUCCUAGGUAUUCGAGCCAAUAAGUUCGG",
    );
    defer query.deinit();

    // reverse complement so we test minus strand search
    query.reverse();
    query.complement();

    {
        var hits = SearchHitList(alphabet.DNA).init(allocator);
        defer hits.deinit();

        var search = try Search(databaseType).init(allocator, &database, .{ .strand = .plus });
        defer search.deinit();

        try search.search(query, &hits);
        try std.testing.expectEqual(@as(usize, 0), hits.list.items.len);
    }
    {
        var hits = SearchHitList(alphabet.DNA).init(allocator);
        defer hits.deinit();

        var search = try Search(databaseType).init(allocator, &database, .{ .strand = .minus });
        defer search.deinit();

        try search.search(query, &hits);
        try std.testing.expectEqual(@as(usize, 1), hits.list.items.len);
    }
    {
        var hits = SearchHitList(alphabet.DNA).init(allocator);
        defer hits.deinit();

        var search = try Search(databaseType).init(allocator, &database, .{ .strand = .both });
        defer search.deinit();

        try search.search(query, &hits);
        try std.testing.expectEqual(@as(usize, 1), hits.list.items.len);
    }
}
