const std = @import("std");
const print = std.debug.print;

const ArgParser = @import("arg_parser.zig").ArgParser;
const Strand = @import("search.zig").Strand;

const Mode = enum {
    dna,
    protein,
};

pub const Args = struct {
    query: [std.fs.MAX_PATH_BYTES:0]u8 = undefined,
    db: [std.fs.MAX_PATH_BYTES:0]u8 = undefined,
    out: [std.fs.MAX_PATH_BYTES:0]u8 = undefined,
    min_identity: f32 = undefined,
    strand: Strand = .both,
    mode: Mode = .dna,

    max_hits: u32 = 1,
    max_rejects: u32 = 16,
};

const ArgsError = error{
    InsufficientArgs,
};

fn printUsage() void {
    print(
        \\Zig-powered nsearch
        \\
        \\Usage:
        \\ nsearchz --query <queryfile> --db <dbfile> --out <outfile> --min-identity <minidentity> [--max-hits <maxaccepts>] [--max-rejects <maxrejects>] [--protein] [--strand <strand>]
        \\
        \\Options:
        \\ --min-identity <minidentity>    Minimum identity threshold (e.g. 0.8).
        \\ --max-hits <maxaccepts>         Maximum number of successful hits reported for one query [default: 1].
        \\ --max-rejects <maxrejects>      Abort after this many candidates were rejected [default: 16].
        \\ --strand <strand>               Strand to search on (plus, minus or both). If minus (or both), queries are reverse complemented [default: both].
        \\ --mode <mode>                   dna or protein [default: dna].
    , .{});
}

pub fn parseArgs(allocator: std.mem.Allocator) !Args {
    const required = [_][]const u8{
        "--query",
        "--db",
        "--out",
        "--min-identity",
    };

    const is_arg_missing = for (required) |req| {
        const found = for (std.os.argv) |arg| {
            if (std.mem.eql(u8, req, std.mem.sliceTo(arg, 0))) break true;
        } else false;

        if (!found) break true;
    } else false;

    if (is_arg_missing) {
        printUsage();
        return ArgsError.InsufficientArgs;
    }

    var args = Args{};
    try ArgParser(Args).parse(allocator, &args);
    return args;
}
