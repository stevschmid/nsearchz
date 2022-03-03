const std = @import("std");

pub const WorkerState = struct {
    const Self = @This();
    pub const State = enum(u8) {
        idle,
        running,
        finished,
    };

    state: std.atomic.Atomic(State) = std.atomic.Atomic(State).init(.idle),

    fn get(self: *const Self) State {
        return self.state.load(.SeqCst);
    }

    fn set(self: *Self, state: State) void {
        _ = self.state.swap(state, .SeqCst);
    }
};

pub fn CreateWorker(comptime Input: type, comptime Output: type, comptime MaxBufferSize: usize) type {
    return struct {
        const Self = @This();

        pub const InputNode = std.TailQueue(Input).Node;
        pub const OutputNode = std.TailQueue(Output).Node;

        const DependencyKind = enum(u8) {
            upstream,
            consumer,
        };

        const Dependency = struct {
            state: *const WorkerState,
            kind: DependencyKind,
        };

        allocator: std.mem.Allocator,
        thread: ?std.Thread = null,

        state: WorkerState = .{},

        input_queue: ?*std.atomic.Queue(Input) = null,
        output_queue: ?*std.atomic.Queue(Output) = null,

        input_buffer: std.TailQueue(Input) = .{},
        output_buffer: std.TailQueue(Output) = .{},

        deps: std.ArrayList(Dependency),

        pub fn init(allocator: std.mem.Allocator) Self {
            return Self{
                .allocator = allocator,
                .deps = std.ArrayList(Dependency).init(allocator),
            };
        }

        pub fn deinit(self: *Self) void {
            self.deps.deinit();
        }

        fn process(self: *Self, comptime HandlerType: type, handler: *HandlerType) !void {
            // we are idle
            self.state.set(.idle);

            // check if we can enter running
            while (true) {
                const satisfied = for (self.deps.items) |dep| {
                    if (dep.kind == .upstream) {
                        if (dep.state.get() != .finished) {
                            break false;
                        }
                    }
                } else true;

                if (satisfied) break;
                std.time.sleep(10 * std.time.ns_per_ms);
            }

            // Enter running mode
            self.state.set(.running);

            try handler.init(self);

            while (true) {
                try handler.loop(self);

                // all consumer satisified?
                const satisfied = for (self.deps.items) |dep| {
                    if (dep.kind == .consumer) {
                        if (dep.state.get() != .finished) {
                            break false;
                        }
                    }
                } else true;

                if (satisfied) break;
                std.time.sleep(10 * std.time.ns_per_ms);
            }

            try handler.deinit(self);

            // flush before we are declare finished
            self.flush();

            // okay we good
            self.state.set(.finished);
        }

        pub fn run(self: *Self, comptime HandlerType: type, handler: *HandlerType) !void {
            const entryFn = struct {
                pub fn fun(_self: *Self, _handler: *HandlerType) !void {
                    try _self.process(HandlerType, _handler);
                }
            }.fun;

            self.thread = try std.Thread.spawn(.{}, entryFn, .{ self, handler });
        }

        pub fn join(self: *Self) void {
            std.debug.assert(self.thread != null);

            self.thread.?.join();
        }

        pub fn addDependency(self: *Self, dependency: Dependency) !void {
            try self.deps.append(dependency);
        }

        pub fn pop(self: *Self) ?*InputNode {
            const node = self.input_buffer.popFirst();
            if (node != null) {
                return node.?;
            }

            // try to fill input queue
            var count: usize = 0;
            while (self.input_queue.?.get()) |in_node| {
                self.input_buffer.append(in_node);
                count += 1;

                if (count >= MaxBufferSize)
                    break;
            }

            // try agane
            return self.input_buffer.popFirst();
        }

        pub fn push(self: *Self, node: *OutputNode) void {
            self.output_buffer.append(node);

            // flush if buffer exceeded
            if (self.output_buffer.len >= MaxBufferSize) {
                self.flush();
            }
        }

        fn flush(self: *Self) void {
            while (self.output_buffer.popFirst()) |out_node| {
                self.output_queue.?.put(out_node);
            }
        }
    };
}
