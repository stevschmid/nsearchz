# A partial reimplementation of [nsearch](https://github.com/stevschmid/nsearch/) in [Zig](https://ziglang.org/)

This is a partial rewrite of nsearch from C++ into Zig, with the primary goal to learn Zig.
Compared to the C++ version, this version offers a limited feature set, but with better performance.

* Search only, no merge or filter capabilities
* FASTA file support (no FASTQ, no gzip support for input)
* ALNOUT file support
* DNA/protein sequences are both supported

It's unlikely that this project will ever meet feature parity with the C++ version of (lib)nsearch.

## Building

You only need Zig. Built in the project root folder with (`-Drelease-fast` is recommended for best performance, omit for debug).

```bash
zig build -Drelease-fast
```

The binary will be located in `./zig-out/bin/nsearchz`

## Usage example

```bash
./nsearchz --query /opt/data/query.fa  --db /opt/data/db.fa --out /tmp/results.alnout --min-identity 0.75 --max-hits 1 --max-rejects 8 --strand plus
```
