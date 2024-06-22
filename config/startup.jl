#!/usr/bin/env -S julia --threads auto --color=yes --compile=all --optimize=3

module Startup

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
