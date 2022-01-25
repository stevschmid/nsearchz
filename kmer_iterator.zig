const std = @import("std");
const alphabet = @import("alphabet.zig");

pub fn KmerIterator(comptime A: type) type {
    return struct {
        const Self = @This();

        const AlphabetInfo = alphabet.AlphabetInfo(A);

        pos: usize,
        letters: []const u8,
        val: AlphabetInfo.KmerType = 0,
        ambiguity_count: usize = 0,

        fn init(letters: []const u8) Self {
            return Self{
                .pos = 0,
                .letters = letters,
            };
        }

        pub fn next(self: *Self) ?AlphabetInfo.KmerType {
            while (self.pos + 1 < AlphabetInfo.NumLettersPerKmer and self.consumeNext()) {
                // advance up to pos - 1 for the initial kmer
            }

            if (!self.consumeNext())
                return null;

            return if (self.ambiguity_count > 0) AlphabetInfo.AmbiguousKmer else self.val;
        }

        fn consumeNext(self: *Self) bool {
            if (self.pos >= self.letters.len) {
                return false;
            }

            // evaluate current letter
            const letterBits = A.mapToBits(self.letters[self.pos]);
            if (letterBits == null) {
                // the next X kmers are considered to be ambiguous
                self.ambiguity_count = AlphabetInfo.NumLettersPerKmer + 1;
            } else {
                // map current letter
                self.val = ((self.val << @bitSizeOf(AlphabetInfo.LetterToBitMapType)) | letterBits.?) & AlphabetInfo.KmerMask;
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

test "basic test" {
    var it = KmerIterator(alphabet.DNA).init("AAAATTTTCCCCGGGG");
    try std.testing.expectEqual(@intCast(alphabet.AlphabetInfo(alphabet.DNA).KmerType, 0b00_00_00_00_10_10_10_10), it.next().?);
    try std.testing.expectEqual(@intCast(alphabet.AlphabetInfo(alphabet.DNA).KmerType, 0b00_00_00_10_10_10_10_01), it.next().?);
    try std.testing.expectEqual(@intCast(alphabet.AlphabetInfo(alphabet.DNA).KmerType, 0b00_00_10_10_10_10_01_01), it.next().?);
    try std.testing.expectEqual(@intCast(alphabet.AlphabetInfo(alphabet.DNA).KmerType, 0b00_10_10_10_10_01_01_01), it.next().?);
    try std.testing.expectEqual(@intCast(alphabet.AlphabetInfo(alphabet.DNA).KmerType, 0b10_10_10_10_01_01_01_01), it.next().?);

    try std.testing.expectEqual(@intCast(alphabet.AlphabetInfo(alphabet.DNA).KmerType, 0b10_10_10_01_01_01_01_11), it.next().?);
    try std.testing.expectEqual(@intCast(alphabet.AlphabetInfo(alphabet.DNA).KmerType, 0b10_10_01_01_01_01_11_11), it.next().?);
    try std.testing.expectEqual(@intCast(alphabet.AlphabetInfo(alphabet.DNA).KmerType, 0b10_01_01_01_01_11_11_11), it.next().?);
    try std.testing.expectEqual(@intCast(alphabet.AlphabetInfo(alphabet.DNA).KmerType, 0b01_01_01_01_11_11_11_11), it.next().?);

    try std.testing.expect(it.next() == null);
}

test "ambiguous nucleotides" {
    var it = KmerIterator(alphabet.DNA).init("AAAATTTTNTNCCCCGGGG");
    // str len = 19, frame = 8, (19-8)+1 = 12 total frames
    try std.testing.expectEqual(@intCast(alphabet.AlphabetInfo(alphabet.DNA).KmerType, 0b0_00_00_00_00_10_10_10_10), it.next().?);
    var i: usize = 0;
    while (i < 10) : (i += 1) {
        try std.testing.expectEqual(@intCast(alphabet.AlphabetInfo(alphabet.DNA).KmerType, 0b1_11_11_11_11_11_11_11_11), it.next().?); // ambiguous kmer
    }
    try std.testing.expectEqual(@intCast(alphabet.AlphabetInfo(alphabet.DNA).KmerType, 0b01_01_01_01_11_11_11_11), it.next().?);
    try std.testing.expect(it.next() == null);
}
