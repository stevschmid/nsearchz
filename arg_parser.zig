const std = @import("std");

// Generic arg parser
pub fn ArgParser(comptime T: type) type {
    return struct {
        pub fn parse(allocator: std.mem.Allocator, values: *T) !void {
            const args = try std.process.argsAlloc(allocator);
            defer std.process.argsFree(allocator, args);

            inline for (@typeInfo(T).Struct.fields) |field| {
                const info = @typeInfo(field.field_type);

                for (args) |arg, index| {
                    if (!std.mem.startsWith(u8, arg, "--"))
                        continue;

                    const arg_field_name = try allocator.dupe(u8, arg[2..]);
                    defer allocator.free(arg_field_name);

                    // flag-name => flag_name
                    std.mem.replaceScalar(u8, arg_field_name, '-', '_');

                    if (!std.mem.eql(u8, field.name, arg_field_name))
                        continue;

                    if (index + 1 >= args.len and info != .Bool) // make sure we have an actual value (if not bool)
                        continue;

                    switch (info) {
                        .Bool => @field(values, field.name) = true,
                        .Int => @field(values, field.name) = try std.fmt.parseInt(field.field_type, args[index + 1], 10),
                        .Float => @field(values, field.name) = try std.fmt.parseFloat(field.field_type, args[index + 1]),
                        .Array => std.mem.copy(u8, &@field(values, field.name), args[index + 1]),
                        .Enum => {
                            inline for (info.Enum.fields) |enum_field| {
                                if (std.mem.eql(u8, enum_field.name, args[index + 1])) {
                                    @field(values, field.name) = @intToEnum(field.field_type, enum_field.value);
                                }
                            }
                        },
                        else => @compileError("Argparse of type " ++ @typeName(
                            field.field_type,
                        ) ++ " not supported"),
                    }
                }
            }
        }
    };
}
