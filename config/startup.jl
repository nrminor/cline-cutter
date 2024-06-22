#!/usr/bin/env -S julia --threads auto --color=yes --compile=all --optimize=3

module Startup

try
    using PrecompileTools
catch err
    @warn "Error importing PrecompileTools: $err. Code loading could not be accelerated. Exiting now."
    exit(0)
end

using PrecompileTools

@setup_workload begin
    @compile_workload begin
        using Pkg
        using CSV
        using DataFrames
        using Logging
        using Pipe
        using PrecompileTools
        using VCFTools
        using VariantCallFormat
    end
end

end
