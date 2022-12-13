#!/usr/bin/env julia

using Pkg
Pkg.activate(@__DIR__)

using ArgParse
using JuliaFormatter

function main(args)
    settings = ArgParseSettings()

    @add_arg_table! settings begin
        "--check", "-c"
        help = "Check to see if formatting is OK"
        action = "store_true"

        "--staged", "-s"
        help = "Format staged files only"
        action = "store_true"
    end

    # apply setting to args
    parsed_args = parse_args(args, settings; as_symbols=true)

    # set overwrite to true; set to false if only checking
    options = Dict(:overwrite => true)
    if parsed_args[:check]
        options[:overwrite] = false
    end

    status = format(pwd(), BlueStyle(); options...) ? 0 : 1
    return exit(status)
end

main(ARGS)
