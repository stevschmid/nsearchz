const std = @import("std");
const assert = std.debug.assert;

pub const DNA = struct {
    pub const NumLettersPerKmer = 8;

    pub fn complement(letter: u8) u8 {
        return switch (letter) {
            'A' => 'T',
            'G' => 'C',
            'C' => 'G',

            'T',
            'U' => 'A',

            'Y' => 'R',
            'R' => 'Y',
            'W' => 'W',
            'S' => 'S',
            'K' => 'M',
            'M' => 'K',

            'D' => 'H',
            'V' => 'B',
            'H' => 'D',
            'B' => 'V',
            'N' => 'N',
            else => letter,
        };
    }

    pub fn mapToBits(letter: u8) ?u2 {
        return switch (letter) {
            'A' => 0b00,
            'C' => 0b01,
            'U',
            'T' => 0b10,
            'G' => 0b11,
            else => null,
        };
    }

    const ScoreMatrixSize = 26;
    const ScoreMatrix = [ScoreMatrixSize][ScoreMatrixSize]i8{
        [_]i8{  2, -4, -4,  2,  0,  0, -4,  2,  0,  0, -4,  0,  2,  2,  0,  0,  0,  2, -4, -4, -4,  2,  2,  0, -4,  0 },
        [_]i8{ -4,  2,  2,  2,  0,  0,  2,  2,  0,  0,  2,  0,  2,  2,  0,  0,  0,  2,  2,  2,  2,  2,  2,  0,  2,  0 },
        [_]i8{ -4,  2,  2, -4,  0,  0, -4,  2,  0,  0, -4,  0,  2,  2,  0,  0,  0, -4,  2, -4, -4,  2, -4,  0,  2,  0 },
        [_]i8{  2,  2, -4,  2,  0,  0,  2,  2,  0,  0,  2,  0,  2,  2,  0,  0,  0,  2,  2,  2,  2,  2,  2,  0,  2,  0 },
        [_]i8{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
        [_]i8{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
        [_]i8{ -4,  2, -4,  2,  0,  0,  2, -4,  0,  0,  2,  0, -4,  2,  0,  0,  0,  2,  2, -4, -4,  2, -4,  0, -4,  0 },
        [_]i8{  2,  2,  2,  2,  0,  0, -4,  2,  0,  0,  2,  0,  2,  2,  0,  0,  0,  2,  2,  2,  2,  2,  2,  0,  2,  0 },
        [_]i8{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
        [_]i8{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
        [_]i8{ -4,  2, -4,  2,  0,  0,  2,  2,  0,  0,  2,  0, -4,  2,  0,  0,  0,  2,  2,  2,  2,  2,  2,  0,  2,  0 },
        [_]i8{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
        [_]i8{  2,  2,  2,  2,  0,  0, -4,  2,  0,  0, -4,  0,  2,  2,  0,  0,  0,  2,  2, -4, -4,  2,  2,  0,  2,  0 },
        [_]i8{  2,  2,  2,  2,  0,  0,  2,  2,  0,  0,  2,  0,  2,  2,  0,  0,  0,  2,  2,  2,  2,  2,  2,  0,  2,  0 },
        [_]i8{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
        [_]i8{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
        [_]i8{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
        [_]i8{  2,  2, -4,  2,  0,  0,  2,  2,  0,  0,  2,  0,  2,  2,  0,  0,  0,  2,  2, -4, -4,  2,  2,  0, -4,  0 },
        [_]i8{ -4,  2,  2,  2,  0,  0,  2,  2,  0,  0,  2,  0,  2,  2,  0,  0,  0,  2,  2, -4, -4,  2, -4,  0,  2,  0 },
        [_]i8{ -4,  2, -4,  2,  0,  0, -4,  2,  0,  0,  2,  0, -4,  2,  0,  0,  0, -4, -4,  2,  2, -4,  2,  0,  2,  0 },
        [_]i8{ -4,  2, -4,  2,  0,  0, -4,  2,  0,  0,  2,  0, -4,  2,  0,  0,  0, -4, -4,  2,  2, -4,  2,  0,  2,  0 },
        [_]i8{  2,  2,  2,  2,  0,  0,  2,  2,  0,  0,  2,  0,  2,  2,  0,  0,  0,  2,  2, -4, -4,  2,  2,  0,  2,  0 },
        [_]i8{  2,  2, -4,  2,  0,  0, -4,  2,  0,  0,  2,  0,  2,  2,  0,  0,  0,  2, -4,  2,  2,  2,  2,  0,  2,  0 },
        [_]i8{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
        [_]i8{ -4,  2,  2,  2,  0,  0, -4,  2,  0,  0,  2,  0,  2,  2,  0,  0,  0, -4,  2,  2,  2,  2,  2,  0,  2,  0 },
        [_]i8{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
    };

    pub fn score(a: u8, b: u8) i8 {
        const idx1: usize = a - 'A';
        const idx2: usize = b - 'A';

        assert(idx1 <= ScoreMatrixSize and idx2 <= ScoreMatrixSize);

        return ScoreMatrix[idx1][idx2];
    }

    pub fn match(a: u8, b: u8) bool {
        return score(a, b) > 0;
    }
};

test "complement" {
    try std.testing.expectEqual(DNA.complement('A'), 'T');
    try std.testing.expectEqual(DNA.complement('G'), 'C');
    try std.testing.expectEqual(DNA.complement('C'), 'G');
    try std.testing.expectEqual(DNA.complement('T'), 'A');
    try std.testing.expectEqual(DNA.complement('U'), 'A');
}

test "bit mapping" {
    try std.testing.expectEqual(DNA.mapToBits('A'), 0b00);

    try std.testing.expectEqual(DNA.mapToBits('C'), 0b01);

    try std.testing.expectEqual(DNA.mapToBits('U'), 0b10);
    try std.testing.expectEqual(DNA.mapToBits('T'), 0b10);

    try std.testing.expectEqual(DNA.mapToBits('G'), 0b11);

    try std.testing.expectEqual(DNA.mapToBits('N'), 0b11);
}

test "scoring" {
    try std.testing.expectEqual(DNA.score('A', 'A'), 2);
    try std.testing.expectEqual(DNA.score('A', 'T'), -4);

    try std.testing.expectEqual(DNA.score('A', 'M'), 2); // M = amino (A, C)
}

test "match" {
    try std.testing.expectEqual(DNA.match('A', 'A'), true);
    try std.testing.expectEqual(DNA.match('C', 'G'), false);
    try std.testing.expectEqual(DNA.match('G', 'G'), true);

    try std.testing.expectEqual(DNA.match('K', 'U'), true); // K = keto (T/U G)
    try std.testing.expectEqual(DNA.match('K', 'M'), false);
}
