#!/usr/bin/env -S julia --threads auto --color=yes --gcthreads 3 --compile=all

module Startup

using PrecompileTools

@setup_workload begin
    @compile_workload begin
        using Pkg
        using CSV
        using DataFrames
        using Logging
        using Pipe
        using VCFTools
        using VariantCallFormat
        using AbbreviatedStackTraces
    end
end

end
