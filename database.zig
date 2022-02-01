const std = @import("std");
const bio = @import("bio/bio.zig");

const Sequence = @import("sequence.zig").Sequence;

pub fn Database(comptime A: type, comptime KmerLength: comptime_int) type {
    return struct {
        const Self = @This();
        pub const kmerInfo = bio.kmer.KmerInfo(A, KmerLength);

        allocator: std.mem.Allocator,
        sequences: []Sequence(A),

        seq_indices: []const usize,
        seq_offset_by_kmer: []const usize,
        seq_count_by_kmer: []const usize,

        kmers: []const kmerInfo.KmerType,
        kmer_offset_by_seq: []const usize,
        kmer_count_by_seq: []const usize,

        // Database takes ownerhip of sequences
        pub fn init(allocator: std.mem.Allocator, sequences: []Sequence(A)) !Self {
            var total_entries: usize = 0;
            var total_unique_entries: usize = 0;

            // to keep track of the unique kmer of a given sequence
            var seq_by_kmer = try allocator.alloc(usize, kmerInfo.MaxKmers);
            std.mem.set(usize, seq_by_kmer, std.math.maxInt(usize)); // set to -1 equivalent
            defer allocator.free(seq_by_kmer);

            var count_by_kmer = try allocator.alloc(usize, kmerInfo.MaxKmers);
            std.mem.set(usize, count_by_kmer, 0);
            defer allocator.free(count_by_kmer);

            for (sequences) |sequence, seq_idx| {
                var kmer_it = bio.kmer.Iterator(kmerInfo).init(sequence.data);
                while (kmer_it.next()) |kmer| {
                    total_entries += 1;

                    // ambiguous?
                    if (kmer == kmerInfo.AmbiguousKmer)
                        continue;

                    // already counted for this sequence?
                    if (seq_by_kmer[kmer] == seq_idx)
                        continue;

                    seq_by_kmer[kmer] = seq_idx;
                    count_by_kmer[kmer] += 1;
                    total_unique_entries += 1;
                }
            }

            // Populate DB
            var kmers = try allocator.alloc(kmerInfo.KmerType, total_entries);
            var kmer_offset_by_seq = try allocator.alloc(usize, sequences.len);
            var kmer_count_by_seq = try allocator.alloc(usize, sequences.len);

            var seq_indices = try allocator.alloc(usize, total_unique_entries);
            var seq_count_by_kmer = try allocator.alloc(usize, kmerInfo.MaxKmers);
            var seq_offset_by_kmer = try allocator.alloc(usize, kmerInfo.MaxKmers);

            // zero counts
            std.mem.set(usize, kmer_count_by_seq, 0);
            std.mem.set(usize, seq_count_by_kmer, 0);

            for (seq_offset_by_kmer) |*seq_offset, kmer| {
                seq_offset.* = if (kmer > 0) seq_offset_by_kmer[kmer - 1] + count_by_kmer[kmer - 1] else 0;
            }

            // Reset tracking for unique kmer within a sequence
            std.mem.set(usize, seq_by_kmer, std.math.maxInt(usize)); // set to -1 equivalent

            var kmer_count: usize = 0;
            for (sequences) |sequence, seq_idx| {
                kmer_offset_by_seq[seq_idx] = kmer_count;

                var kmer_it = bio.kmer.Iterator(kmerInfo).init(sequence.data);
                while (kmer_it.next()) |kmer| {
                    kmers[kmer_count] = kmer;
                    kmer_count += 1;

                    // ambiguous?
                    if (kmer == kmerInfo.AmbiguousKmer)
                        continue;

                    // already counted for this sequence?
                    if (seq_by_kmer[kmer] == seq_idx)
                        continue;

                    seq_by_kmer[kmer] = seq_idx;

                    seq_indices[seq_offset_by_kmer[kmer] + seq_count_by_kmer[kmer]] = seq_idx;
                    seq_count_by_kmer[kmer] += 1;
                }

                kmer_count_by_seq[seq_idx] = kmer_count - kmer_offset_by_seq[seq_idx];
            }

            return Self{
                .allocator = allocator,
                .sequences = sequences,

                .seq_indices = seq_indices,
                .seq_offset_by_kmer = seq_offset_by_kmer,
                .seq_count_by_kmer = seq_count_by_kmer,

                .kmers = kmers,
                .kmer_offset_by_seq = kmer_offset_by_seq,
                .kmer_count_by_seq = kmer_count_by_seq,
            };
        }

        pub fn deinit(self: *Self) void {
            self.allocator.free(self.seq_indices);
            self.allocator.free(self.seq_offset_by_kmer);
            self.allocator.free(self.seq_count_by_kmer);

            self.allocator.free(self.kmers);
            self.allocator.free(self.kmer_offset_by_seq);
            self.allocator.free(self.kmer_count_by_seq);

            for (self.sequences) |*sequence| sequence.deinit();
            self.allocator.free(self.sequences);
        }
    };
}
