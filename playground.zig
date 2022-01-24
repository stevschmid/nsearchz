const std = @import("std");
const print = std.debug.print;

// pub const Alphabet = enum {
//     DNA,
//     Protein,

//     pub fn numBits(comptime self: Alphabet) u8 {
//         return switch (self) {
//             Alphabet.DNA => 2,
//             Alphabet.Protein => 4,
//         };
//     }

//     pub fn complement(comptime self: Alphabet, letter: u8) u8 {
//         return switch (self) {
//             Alphabet.DNA => {
//                 return switch (letter) {
//                     'A' => 'T',
//                     else => letter,
//                 };
//             },
//             else => unreachable,
//         };
//     }
// };
// const Alphabet = enum { DNA, Protein };

// pub fn complement(comptime A: Alphabet, letter: u8) u8 {
//     return switch(A) {
//         .DNA => letter,
//         .Protein => letter,
//     };
// }

const DNA = struct {
    pub const num_bits = 2;

    pub fn bitMap(letter: u8) u8 {
        _ = letter;
        return 1;
    }
};

const Protein = struct {
    pub const num_bits = 4;

    pub fn bitMap(letter: u8) u8 {
        _ = letter;
        return 2;
    }
};

pub fn doThings(comptime A: type) void {
    print("Bitmap {d} num_bits: {d}\n", .{A.bitMap('A'), A.num_bits});
}

pub fn main() !void {
    // const iter = Iterator(u8);
    // iter.nextFn();
    // _ = iter;

    // var dna: Alphabet = Alphabet.DNA;
    // _ = dna;

    // print("Complement {d}\n", .{DNA.num_bits});

    // const dna: DNA = DNA {};

    // print("Complement {c}\n", .{complement(Alphabet.DNA, 'T')});
    // print("Bitmap {d}\n", .{DNA.bitMap('T')});
    doThings(DNA);
    doThings(struct {});
}
