const std = @import("std");
const alphabet = @import("alphabet.zig");

pub fn KmerInfo(comptime A: type, comptime NumLettersPerKmer: comptime_int) type {
    return struct {
        pub const NumLettersPerKmer = NumLettersPerKmer;
        pub const Alphabet = A;
        pub const LetterToBitMapType = @typeInfo(@typeInfo(@TypeOf(A.mapToBits)).Fn.return_type.?).Optional.child; // ?u2 -> u2
        pub const NumBitsPerLetter = @bitSizeOf(LetterToBitMapType);
        pub const NumBitsPerKmer = (NumLettersPerKmer * NumBitsPerLetter) + 1; // +1 for ambiguous kmer
        pub const KmerType = @Type(.{
            .Int = .{
                .signedness = .unsigned,
                .bits = NumBitsPerKmer, // save ambiguous kmer
            },
        });
        pub const MaxKmers = (1 << (NumBitsPerKmer - 1)) + 1;
        pub const AmbiguousKmer = (1 << NumBitsPerKmer) - 1;
        pub const KmerMask = AmbiguousKmer >> 1;
    };
}

pub fn Iterator(comptime kmerInfo: type) type {
    return struct {
        const Self = @This();
        const A = kmerInfo.Alphabet;

        pos: usize,
        letters: []const u8,
        val: kmerInfo.KmerType = 0,
        ambiguity_count: usize = 0,

        pub fn init(letters: []const u8) Self {
            return Self{
                .pos = 0,
                .letters = letters,
            };
        }

        pub fn num_total(self: *Self) usize {
            if (kmerInfo.NumLettersPerKmer > self.letters.len) {
                return 0;
            } else {
                return self.letters.len - kmerInfo.NumLettersPerKmer + 1;
            }
        }

        pub fn next(self: *Self) ?kmerInfo.KmerType {
            while (self.pos + 1 < kmerInfo.NumLettersPerKmer and self.consumeNext()) {
                // advance up to pos - 1 for the initial kmer
            }

            if (!self.consumeNext())
                return null;

            return if (self.ambiguity_count > 0) kmerInfo.AmbiguousKmer else self.val;
        }

        fn consumeNext(self: *Self) bool {
            if (self.pos >= self.letters.len) {
                return false;
            }

            // evaluate current letter
            const letterBits = A.mapToBits(self.letters[self.pos]);
            if (letterBits == null) {
                // the next X kmers are considered to be ambiguous
                self.ambiguity_count = kmerInfo.NumLettersPerKmer + 1;
            } else {
                // map current letter
                self.val = ((self.val << @bitSizeOf(kmerInfo.LetterToBitMapType)) | letterBits.?) & kmerInfo.KmerMask;
            }

            // advance
            self.pos += 1;

            // in ambiguous region?
            if (self.ambiguity_count > 0) {
                self.ambiguity_count -= 1;
            }

            return true;
        }
    };
}

test "info" {
    const a = KmerInfo(alphabet.DNA, 8);

    try std.testing.expectEqual(8, a.NumLettersPerKmer);
    try std.testing.expectEqual(17, a.NumBitsPerKmer); //2*8+1
    try std.testing.expectEqual(u2, a.LetterToBitMapType);
    try std.testing.expectEqual(u17, a.KmerType);
    try std.testing.expectEqual(65537, a.MaxKmers);
    try std.testing.expectEqual(0b1_11_11_11_11_11_11_11_11, a.AmbiguousKmer); // ambiguous flag active
    try std.testing.expectEqual(0b0_11_11_11_11_11_11_11_11, a.KmerMask);
}

test "basic test" {
    const kmerInfo = KmerInfo(alphabet.DNA, 3);

    var it = Iterator(kmerInfo).init("ATCGGG");
    try std.testing.expectEqual(@as(kmerInfo.KmerType, 0b00_10_01), it.next().?);
    try std.testing.expectEqual(@as(kmerInfo.KmerType, 0b10_01_11), it.next().?);
    try std.testing.expectEqual(@as(kmerInfo.KmerType, 0b01_11_11), it.next().?);
    try std.testing.expectEqual(@as(kmerInfo.KmerType, 0b11_11_11), it.next().?);
    try std.testing.expect(it.next() == null);

    try std.testing.expectEqual(@as(usize, 4), it.num_total());
}

test "ambiguous nucleotides" {
    const kmerInfo = KmerInfo(alphabet.DNA, 3);
    var it = Iterator(kmerInfo).init("ATCGNGTTNAAGN");

    try std.testing.expectEqual(@as(kmerInfo.KmerType, 0b00_10_01), it.next().?); // ATC
    try std.testing.expectEqual(@as(kmerInfo.KmerType, 0b10_01_11), it.next().?); // TCG

    try std.testing.expectEqual(@as(kmerInfo.KmerType, 0b1_11_11_11), it.next().?); // CGN ambiguous
    try std.testing.expectEqual(@as(kmerInfo.KmerType, 0b1_11_11_11), it.next().?); // GNG skipped
    try std.testing.expectEqual(@as(kmerInfo.KmerType, 0b1_11_11_11), it.next().?); // NGT skipped

    try std.testing.expectEqual(@as(kmerInfo.KmerType, 0b11_10_10), it.next().?); // GTT

    try std.testing.expectEqual(@as(kmerInfo.KmerType, 0b1_11_11_11), it.next().?); // TTN skipped
    try std.testing.expectEqual(@as(kmerInfo.KmerType, 0b1_11_11_11), it.next().?); // TNA skipped
    try std.testing.expectEqual(@as(kmerInfo.KmerType, 0b1_11_11_11), it.next().?); // NAA skipped

    try std.testing.expectEqual(@as(kmerInfo.KmerType, 0b00_00_11), it.next().?); // AAG

    try std.testing.expectEqual(@as(kmerInfo.KmerType, 0b1_11_11_11), it.next().?); // AGN skipped

    try std.testing.expect(it.next() == null);
}

test "too short" {
    const kmerInfo = KmerInfo(alphabet.DNA, 4);
    var it = Iterator(kmerInfo).init("ATT");

    try std.testing.expectEqual(@as(usize, 0), it.num_total());
    try std.testing.expect(it.next() == null);
}

test "just right" {
    const kmerInfo = KmerInfo(alphabet.DNA, 4);
    var it = Iterator(kmerInfo).init("ATTG");

    try std.testing.expectEqual(@as(usize, 1), it.num_total());
    try std.testing.expectEqual(@as(kmerInfo.KmerType, 0b00_10_10_11), it.next().?);
    try std.testing.expect(it.next() == null);
}
