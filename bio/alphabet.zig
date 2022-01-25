pub const DNA = @import("alphabet/dna.zig").DNA;
pub const Protein = @import("alphabet/protein.zig").Protein;

pub fn AlphabetInfo(comptime A: type) type {
    return struct{
        pub const NumLettersPerKmer = A.NumLettersPerKmer;

        pub const LetterToBitMapType = @typeInfo(@typeInfo(@TypeOf(A.mapToBits)).Fn.return_type.?).Optional.child; // ?u2 -> u2
        pub const NumBitsPerLetter = @bitSizeOf(LetterToBitMapType);
        pub const NumBitsPerKmer = (NumLettersPerKmer * NumBitsPerLetter) + 1; // +1 for ambiguous kmer
        pub const KmerType = @Type(.{ .Int = .{
            .signedness = .unsigned,
            .bits = NumBitsPerKmer, // save ambiguous kmer
        } });
        pub const MaxKmers = 1 << NumBitsPerKmer;
        pub const AmbiguousKmer = MaxKmers - 1;
        pub const KmerMask = AmbiguousKmer >> 1; // remove ambiguous bit
    };
}

test "check" {
    const a = AlphabetInfo(DNA);
    const std = @import("std");

    try std.testing.expectEqual(8, a.NumLettersPerKmer);
    try std.testing.expectEqual(17, a.NumBitsPerKmer); //2*8+1
    try std.testing.expectEqual(u2, a.LetterToBitMapType);
    try std.testing.expectEqual(u17, a.KmerType);
    try std.testing.expectEqual(131_072, a.MaxKmers);
    try std.testing.expectEqual(0b1_11_11_11_11_11_11_11_11, a.AmbiguousKmer); // ambiguous flag active
    try std.testing.expectEqual(0b0_11_11_11_11_11_11_11_11, a.KmerMask);
}
