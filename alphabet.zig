pub const DNA = @import("alphabet/dna.zig").DNA;
pub const Protein = @import("alphabet/protein.zig").Protein;

pub fn AlphabetInfo(comptime A: type) type {
    return struct{
        pub const NumLettersPerKmer = A.NumLettersPerKmer;

        pub const LetterToBitMapType = @typeInfo(@typeInfo(@TypeOf(A.mapToBits)).Fn.return_type.?).Optional.child; // ?u2 -> u2
        pub const NumBitsPerLetter = @bitSizeOf(LetterToBitMapType);
        pub const NumBitsPerKmer = NumLettersPerKmer * NumBitsPerLetter;
        pub const KmerType = @Type(.{ .Int = .{
            .signedness = .unsigned,
            .bits = NumBitsPerKmer,
        } });
        pub const KmerMask =  ~@intCast(KmerType, 0);
        pub const MaxKmers = 1 << NumBitsPerKmer;
    };
}
