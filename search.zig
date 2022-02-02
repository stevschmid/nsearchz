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

pub const SearchOptions = struct {
    max_accepts: u32 = 1,
    max_rejects: u32 = 16,
    min_identity: f32 = 0.75,
};

pub const SearchHit = struct {
    const Self = @This();

    allocator: std.mem.Allocator,
    db_seq_id: usize,
    cigar: Cigar,

    pub fn init(allocator: std.mem.Allocator) Self {
        return Self{
            .allocator = allocator,
            .cigar = Cigar.init(allocator),
            .db_seq_id = undefined,
        };
    }

    pub fn deinit(self: *Self) void {
        self.cigar.deinit();
    }
};

pub const SearchHitList = utils.ArrayListDeinitWrapper(SearchHit);
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

        pub fn search(self: *Self, query: Sequence(A), hits: *SearchHitList) !void {
            const min_hsp_length = std.math.min(DefaultMinHspLength, query.data.len / 2);
            const max_hsp_join_distance = DefaultMaxHSPJoinDistance;
            _ = min_hsp_length;
            _ = max_hsp_join_distance;

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
                        _ = try self.banded_align.process(query, seq, .forward, part.hsp.end_one + 1, part.hsp.end_two + 1, next_part.hsp.start_one, next_part.hsp.start_two, &cigar);
                        try final_cigar.appendOther(cigar);
                    }

                    // Align last HSP's end to whole sequences end
                    try final_cigar.appendOther(last_part.cigar);
                    _ = try self.banded_align.process(query, seq, .forward, last_part.hsp.end_one + 1, last_part.hsp.end_two + 1, null, null, &cigar);
                    try final_cigar.appendOther(cigar);

                    var accept = (final_cigar.identity() >= self.options.min_identity);
                    if (accept) {
                        var hit = SearchHit.init(hits.allocator);
                        hit.db_seq_id = seq_id;
                        try hit.cigar.appendOther(final_cigar);

                        try hits.list.append(hit);

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

    var database = try databaseType.init(allocator, sequences.list.toOwnedSlice());
    defer database.deinit();

    {
        var search = try Search(databaseType).init(allocator, &database, .{});
        defer search.deinit();

        var hits = SearchHitList.init(allocator);
        defer hits.deinit();

        try search.search(query, &hits);

        try std.testing.expectEqual(@as(usize, 1), hits.list.items.len);
        try std.testing.expectEqualStrings("DB1", database.sequences[hits.list.items[0].db_seq_id].identifier);
    }

    // try accepts 2, other sequence is still low
    {
        var search = try Search(databaseType).init(allocator, &database, .{ .max_accepts = 2 });
        defer search.deinit();

        var hits = SearchHitList.init(allocator);
        defer hits.deinit();

        try search.search(query, &hits);

        // still 1
        try std.testing.expectEqual(@as(usize, 1), hits.list.items.len);
    }

    // accept two, but lower identity threshold
    {
        var search = try Search(databaseType).init(allocator, &database, .{ .max_accepts = 2, .min_identity = 0.5 });
        defer search.deinit();

        var hits = SearchHitList.init(allocator);
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

    var database = try databaseType.init(allocator, sequences.list.toOwnedSlice());
    defer database.deinit();

    var search = try Search(databaseType).init(allocator, &database, .{});
    defer search.deinit();

    var query1 = try Sequence(alphabet.DNA).init(allocator, "Query 1 ", "TGAGACGATGCAAA");
    defer query1.deinit();

    var hits1 = SearchHitList.init(allocator);
    defer hits1.deinit();

    try search.search(query1, &hits1);
    try std.testing.expectEqual(@as(usize, 1), hits1.list.items.len);

    var query2 = try Sequence(alphabet.DNA).init(allocator, "Query 1 ", "CGTTATATTCGGAGACCTAT");
    defer query2.deinit();

    var hits2 = SearchHitList.init(allocator);
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

    var database = try databaseType.init(allocator, sequences.list.toOwnedSlice());
    defer database.deinit();

    var search = try Search(databaseType).init(allocator, &database, .{ .min_identity = 0.5 });
    defer search.deinit();

    var hits = SearchHitList.init(allocator);
    defer hits.deinit();

    try search.search(query, &hits);
    try std.testing.expectEqual(@as(usize, 1), hits.list.items.len);

    const hit = hits.list.items[0];

    // expect nsearch output
    try std.testing.expectApproxEqAbs(@as(f32, 0.52), hit.cigar.identity(), 0.01);
    try std.testing.expectEqualStrings("92D4X1=4X2=1X1=1X1=4X3=2X4=2X1=1X1=4X1=2X1=7X3=1X1=1X1=1X1=2X2=1X1=1X1=1X3=2X8=1X1=1D3=19X2=7X2=1X1=9X2=2X1=1X2=2X2=4X2=3X1=3X3=2X1=2X2=1X1=1X1=5X2=2X1=1X1=2X2=2D3=1D1=1D1=1D2X1=1D5X1=2X1=1D2=2X2=2X1=9X2=2X1=2X1=4X3=5X1=1X1=2X1=3X2=1X2=1X1=1X3=4X1=1X4=1D2X1=5X2=3X1=4X1=5X3=3X1=1X3=2X2=6X1=2X1=1D1=5D3=1X1=2X1=2D1X3=3X1=5X1=1X1=2X4=2X1=3X1=3X1=3X3=5X1=3X3=4X4=2X1=7X2=9X1=1X4=2X1=2X1=1X2=1X3=1X2=1X1=2X1=1X1=4X2=1X1=1X3=1D3=2X3=2X1=2X1=2I2X2=1X1=3I1X1=2X2=1X1=2X1=3X3=3X1=2X1=1X1=4X2=1I2X1=2X1=1X2=1X2=1X1=1X1=1I3X1=3I1=2X1=3X1=1X1=1I2X3=2X1=3I3=3X3=2X2=1X1=1X1=2X3=3X1=1X1=5X2=1X1=2X1=1X1=2X1=1X1=2X3=2X1=2X1=2X2=9X1=5X1=2X1=1X1=1X2=2X4=1D3=1X1=6D2=1X1=1X3=1X1=3D1=4D3=1X1=1X2=1X5=2X3=1X14=1X2=1X4=1X13=2I1=2X5=2I2=2X17=1X9=1X2=2X12=1X1=1D2=1X6D1=1X3=1X10=1X3=1X7=1X20=1X3=1X21=1X2=3X2=1X2=1X1=1X7=1X25=1X39=1X2=1X11=1X4=2X8=2X5=1X2=4D2=1X1=1X6=2X8=1X1=1X7=1X1=1X1=2X7=7D2=1X2=1X9=1X1D7=1X11=1X10=1X3=1X5=2X9=1X18=1X20=1X17=2X5=1X2=1X12=1X9=1X5=1X4=2X19=1I2=3I2=1X1=1X5=1X1=2X2=2X1=3X1=1X3=8X2=2X1D2=2X1=3X2=3X1=3X2=2X2=8X1=2X3=1X2=2X2=3X1=2X1=2X3I3=5X2=2X2=3X2=1X1=1X1=1X2=6X3=2X1=2X1=7X2=7X1=1X1=1X2=1X1=1X2=1X1=2X4=8I1=3I1=1X1=2I2=1X1=1X1=6X1=1X2=2X1=1X1=5X1=1X1=1X1=3X1=3X1=4X1=1X1=2X5I1=1X2=1X1=1X2=2X2=1X1=2X1=1X1=1X1=2X2=1X1=2X1=3X1=1X2=4X1=1X1=3X3=1X3=3X1=11X2=1X1=2X2=2X2=1X1=3D1=1D1=5X1=1X1=3X1=2X2=1X1=3X3=1X1=5X1=4X2=4X2=5X1=2X1=1X1=2X1=2X2=6X1=1I3=6X1=2X1=3X1=3X1=2X2=5X1D4=5X3=1X4=1X1=1X2=3X1=2X2=1X1=2X1=1X1=2X2=2X2=3X1=1X2=3X1=11X1=2X2=3X2=1X1=28D", hit.cigar.str());
}
