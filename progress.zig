const std = @import("std");
const utils = @import("utils.zig");

const CountConverter = utils.Converter{
    .units = &[_][]const u8{ "", "k", "M", "G" },
    .dividers = &[_]usize{ 1, 1_000, 1_000_000, 1_000_000_000 },
};

const ByteConverter = utils.Converter{
    .units = &[_][]const u8{ "B", "kB", "MB", "GB" },
    .dividers = &[_]usize{ 1, 1024, 1024 * 1024, 1024 * 1024 * 1024 },
};

pub fn Progress(comptime Stages: type) type {
    return struct {
        const Unit = enum {
            counts,
            bytes,
        };

        const Stage = struct {
            label: []const u8,
            unit: Unit,
            value: usize = 0,
            max: usize = 0,
        };

        const Self = @This();

        const Indexer = std.enums.EnumIndexer(Stages);
        stages: [Indexer.count]Stage = undefined,
        active_index: usize = 0,
        out: std.fs.File = std.io.getStdOut(),

        timer: std.time.Timer = undefined,

        prev_stage: ?*Stage = undefined,
        prev_print_timestamp: u64 = undefined,

        pub fn start(self: *Self) !void {
            self.timer = try std.time.Timer.start();
            self.prev_stage = null;
            self.prev_print_timestamp = 0;
        }

        pub fn finish(self: *Self) void {
            self.write("\n", .{});
        }

        pub fn add(self: *Self, stage: Stages, label: []const u8, unit: Unit) void {
            self.stages[Indexer.indexOf(stage)] = .{
                .label = label,
                .unit = unit,
            };
        }

        pub fn activate(self: *Self, stage: Stages) void {
            const new_index = Indexer.indexOf(stage);
            if (new_index != self.active_index)
                self.write("\n", .{});

            self.active_index = new_index;
        }

        pub fn set(self: *Self, stage: Stages, value: usize, max: usize) void {
            const index = Indexer.indexOf(stage);
            const st = &self.stages[index];
            st.value = value;
            st.max = max;

            if (self.active_index == index) {
                self.print(st);
            }
        }

        fn print(self: *Self, stage: *Stage) void {
            const now = self.timer.read();

            var refresh: bool = false;
            refresh = refresh or (self.prev_print_timestamp + 50 * std.time.ns_per_ms < now);
            refresh = refresh or (self.prev_stage == null or self.prev_stage.? != stage);
            refresh = refresh or (stage.value >= stage.max);

            if (!refresh) {
                return;
            }

            self.prev_print_timestamp = now;
            self.prev_stage = stage;

            const ratio = if (stage.max == 0) 1.0 else (@intToFloat(f32, stage.value) / @intToFloat(f32, stage.max));
            const percent = ratio * 100.0;

            const result = switch (stage.unit) {
                .counts => CountConverter.convert(stage.value),
                .bytes => ByteConverter.convert(stage.value),
            };

            self.write("{s: <20}: {d:3.0}% ({d}{s})" ++ (" " ** 40) ++ "\r", .{ stage.label, percent, result.value, result.unit });
        }

        fn write(self: *Self, comptime fmt: anytype, args: anytype) void {
            self.out.writer().print(fmt, args) catch std.process.exit(1);
        }
    };
}
