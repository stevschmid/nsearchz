const std = @import("std");
const print = std.debug.print;
const ArrayList = std.ArrayList;
const Allocator = std.mem.Allocator;

const utils = @import("utils.zig");

const Sequence = @import("sequence.zig").Sequence;
const FastaReader = @import("fasta_reader.zig").FastaReader;

const bio = @import("bio/bio.zig");
const alphabet = bio.alphabet;

const alphabetChosen = alphabet.DNA;
const kmerInfo = bio.kmer.KmerInfo(alphabetChosen, 8);

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

    // counts by kmer
    var count_by_kmer = try allocator.alloc(usize, kmerInfo.MaxKmers);
    std.mem.set(usize, count_by_kmer, 0);
    defer allocator.free(count_by_kmer);

    // to keep track of the unique kmer of a given sequence
    var seq_by_kmer = try allocator.alloc(usize, kmerInfo.MaxKmers);
    std.mem.set(usize, seq_by_kmer, std.math.maxInt(usize)); // set to -1 equivalent
    defer allocator.free(seq_by_kmer);

    var total_entries: usize = 0;
    var total_unique_entries: usize = 0;

    const sequences = reader.sequences;

    bench_start = std.time.milliTimestamp();
    for (sequences.items) |sequence, seq_idx| {
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

    // Calculate indices
    var seq_offset_by_kmer = try allocator.alloc(usize, kmerInfo.MaxKmers);
    defer allocator.free(seq_offset_by_kmer);

    for (seq_offset_by_kmer) |*seq_offset, kmer| {
        seq_offset.* = if (kmer > 0) seq_offset_by_kmer[kmer - 1] + count_by_kmer[kmer - 1] else 0;
    }

    // Reset tracking for unique kmer within a sequence
    std.mem.set(usize, seq_by_kmer, std.math.maxInt(usize)); // set to -1 equivalent

    // Populate DB
    var kmer_offset_by_seq = try allocator.alloc(usize, sequences.items.len);
    defer allocator.free(kmer_offset_by_seq);

    var kmer_count_by_seq = try allocator.alloc(usize, sequences.items.len);
    defer allocator.free(kmer_count_by_seq);

    var kmers = try allocator.alloc(kmerInfo.KmerType, total_entries);
    defer allocator.free(kmers);

    var seq_indices = try allocator.alloc(usize, total_unique_entries);
    defer allocator.free(seq_indices);

    var seq_count_by_kmer = try allocator.alloc(usize, kmerInfo.MaxKmers);
    std.mem.set(usize, seq_count_by_kmer, 0);
    defer allocator.free(seq_count_by_kmer);

    var kmer_count: usize = 0;
    for (sequences.items) |sequence, seq_idx| {
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

    print("Indexing took {}ms\n", .{std.time.milliTimestamp() - bench_start});

    print("Total Entries {} Unique {}\n", .{total_entries, total_unique_entries});
    print("{b:0>17}\n", .{kmers[kmer_offset_by_seq[1]]});

    // for (count_by_kmer) |count, kmer| {
    //     if (count > 0) {
    //         print("Kmer {b:0>16}: {}\n", .{ kmer, count });
    //     }
    // }
}
