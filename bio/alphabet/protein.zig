const std = @import("std");
const assert = std.debug.assert;

pub const Protein = struct {
    pub fn complement(letter: u8) u8 {
        return letter;
    }

    const BitMapping = [_]?u4{
        0b0000, // 'A'
        null,   // 'B' ambiguous/invalid
        0b0011, // 'C'
        0b0100, // 'D'
        0b0100, // 'E'
        0b1111, // 'F'
        0b0101, // 'G'
        0b0110, // 'H'
        0b0111, // 'I'
        null,   // 'J' ambiguous/invalid
        0b1001, // 'K'
        0b1000, // 'L'
        0b1010, // 'M'
        0b0010, // 'N'
        null,   // 'O' ambiguous/invalid
        0b1011, // 'P'
        0b0100, // 'Q'
        0b0001, // 'R'
        0b1100, // 'S'
        0b1101, // 'T'
        null,   // 'U' ambiguous/invalid
        0b0111, // 'V'
        0b1110, // 'W'
        null,   // 'X' ambiguous/invalid
        0b1111, // 'Y'
        null,   // 'Z' ambiguous/invalid
    };

    pub fn mapToBits(aa: u8) ?u4 {
        const idx: usize = aa - 'A';
        assert(idx <= BitMapping.len);

        return BitMapping[idx];
    }

    // BLOSUM62
    const ScoreMatrixSize = 26;
    const ScoreMatrix = [ScoreMatrixSize][ScoreMatrixSize]i8{
        // A,   B,   C,   D,   E,   F,   G,   H,   I,   J,   K,   L,   M,   N,   O,   P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z
        [_]i8{  4,  -2,   0,  -2,  -1,  -2,   0,  -2,  -1,   0,  -1,  -1,  -1,  -2,   0,  -1,  -1,  -1,   1,   0,   0,   0,  -3,   0,  -2,  -1}, // A
        [_]i8{ -2,   4,  -3,   4,   1,  -3,  -1,   0,  -3,   0,   0,  -4,  -3,   3,   0,  -2,   0,  -1,   0,  -1,   0,  -3,  -4,  -1,  -3,   1}, // B
        [_]i8{  0,  -3,   9,  -3,  -4,  -2,  -3,  -3,  -1,   0,  -3,  -1,  -1,  -3,   0,  -3,  -3,  -3,  -1,  -1,   0,  -1,  -2,  -2,  -2,  -3}, // C
        [_]i8{ -2,   4,  -3,   6,   2,  -3,  -1,  -1,  -3,   0,  -1,  -4,  -3,   1,   0,  -1,   0,  -2,   0,  -1,   0,  -3,  -4,  -1,  -3,   1}, // D
        [_]i8{ -1,   1,  -4,   2,   5,  -3,  -2,   0,  -3,   0,   1,  -3,  -2,   0,   0,  -1,   2,   0,   0,  -1,   0,  -2,  -3,  -1,  -2,   4}, // E
        [_]i8{ -2,  -3,  -2,  -3,  -3,   6,  -3,  -1,   0,   0,  -3,   0,   0,  -3,   0,  -4,  -3,  -3,  -2,  -2,   0,  -1,   1,  -1,   3,  -3}, // F
        [_]i8{  0,  -1,  -3,  -1,  -2,  -3,   6,  -2,  -4,   0,  -2,  -4,  -3,   0,   0,  -2,  -2,  -2,   0,  -2,   0,  -3,  -2,  -1,  -3,  -2}, // G
        [_]i8{ -2,   0,  -3,  -1,   0,  -1,  -2,   8,  -3,   0,  -1,  -3,  -2,   1,   0,  -2,   0,   0,  -1,  -2,   0,  -3,  -2,  -1,   2,   0}, // H
        [_]i8{ -1,  -3,  -1,  -3,  -3,   0,  -4,  -3,   4,   0,  -3,   2,   1,  -3,   0,  -3,  -3,  -3,  -2,  -1,   0,   3,  -3,  -1,  -1,  -3}, // I
        [_]i8{  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0}, // J
        [_]i8{ -1,   0,  -3,  -1,   1,  -3,  -2,  -1,  -3,   0,   5,  -2,  -1,   0,   0,  -1,   1,   2,   0,  -1,   0,  -2,  -3,  -1,  -2,   1}, // K
        [_]i8{ -1,  -4,  -1,  -4,  -3,   0,  -4,  -3,   2,   0,  -2,   4,   2,  -3,   0,  -3,  -2,  -2,  -2,  -1,   0,   1,  -2,  -1,  -1,  -3}, // L
        [_]i8{ -1,  -3,  -1,  -3,  -2,   0,  -3,  -2,   1,   0,  -1,   2,   5,  -2,   0,  -2,   0,  -1,  -1,  -1,   0,   1,  -1,  -1,  -1,  -1}, // M
        [_]i8{ -2,   3,  -3,   1,   0,  -3,   0,   1,  -3,   0,   0,  -3,  -2,   6,   0,  -2,   0,   0,   1,   0,   0,  -3,  -4,  -1,  -2,   0}, // N
        [_]i8{  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0}, // O
        [_]i8{ -1,  -2,  -3,  -1,  -1,  -4,  -2,  -2,  -3,   0,  -1,  -3,  -2,  -2,   0,   7,  -1,  -2,  -1,  -1,   0,  -2,  -4,  -2,  -3,  -1}, // P
        [_]i8{ -1,   0,  -3,   0,   2,  -3,  -2,   0,  -3,   0,   1,  -2,   0,   0,   0,  -1,   5,   1,   0,  -1,   0,  -2,  -2,  -1,  -1,   3}, // Q
        [_]i8{ -1,  -1,  -3,  -2,   0,  -3,  -2,   0,  -3,   0,   2,  -2,  -1,   0,   0,  -2,   1,   5,  -1,  -1,   0,  -3,  -3,  -1,  -2,   0}, // R
        [_]i8{  1,   0,  -1,   0,   0,  -2,   0,  -1,  -2,   0,   0,  -2,  -1,   1,   0,  -1,   0,  -1,   4,   1,   0,  -2,  -3,   0,  -2,   0}, // S
        [_]i8{  0,  -1,  -1,  -1,  -1,  -2,  -2,  -2,  -1,   0,  -1,  -1,  -1,   0,   0,  -1,  -1,  -1,   1,   5,   0,   0,  -2,   0,  -2,  -1}, // T
        [_]i8{  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0}, // U
        [_]i8{  0,  -3,  -1,  -3,  -2,  -1,  -3,  -3,   3,   0,  -2,   1,   1,  -3,   0,  -2,  -2,  -3,  -2,   0,   0,   4,  -3,  -1,  -1,  -2}, // V
        [_]i8{ -3,  -4,  -2,  -4,  -3,   1,  -2,  -2,  -3,   0,  -3,  -2,  -1,  -4,   0,  -4,  -2,  -3,  -3,  -2,   0,  -3,  11,  -2,   2,  -3}, // W
        [_]i8{  0,  -1,  -2,  -1,  -1,  -1,  -1,  -1,  -1,   0,  -1,  -1,  -1,  -1,   0,  -2,  -1,  -1,   0,   0,   0,  -1,  -2,  -1,  -1,  -1}, // X
        [_]i8{ -2,  -3,  -2,  -3,  -2,   3,  -3,   2,  -1,   0,  -2,  -1,  -1,  -2,   0,  -3,  -1,  -2,  -2,  -2,   0,  -1,   2,  -1,   7,  -2}, // Y
        [_]i8{ -1,   1,  -3,   1,   4,  -3,  -2,   0,  -3,   0,   1,  -3,  -1,   0,   0,  -1,   3,   0,   0,  -1,   0,  -2,  -3,  -1,  -2,   4}, // Z
    };

    pub fn score(a: u8, b: u8) i8 {
        var idx1: usize = a - 'A';
        var idx2: usize = b - 'A';

        assert(idx1 <= ScoreMatrixSize and idx2 <= ScoreMatrixSize);

        return ScoreMatrix[idx1][idx2];
    }

    pub fn match(a: u8, b: u8) bool {
        return score(a, b) >= 4;
    }
};

test "complement" {
    try std.testing.expectEqual(Protein.complement('A'), 'A');
}

test "bit mapping" {
    try std.testing.expectEqual(Protein.mapToBits('B'), null);
    try std.testing.expectEqual(Protein.mapToBits('C'), 0b0011);
    try std.testing.expectEqual(Protein.mapToBits('Y'), 0b1111);
}

test "scoring" {
    try std.testing.expectEqual(Protein.score('Q', 'S'), 0);
    try std.testing.expectEqual(Protein.score('L', 'P'), -3);
    try std.testing.expectEqual(Protein.score('T', 'T'), 5);
}

test "match" {
    try std.testing.expectEqual(Protein.match('L', 'L'), true);
    try std.testing.expectEqual(Protein.match('I', 'V'), false);
    try std.testing.expectEqual(Protein.match('I', 'I'), true);
}
