#!/usr/bin/env -S julia --threads auto --history-file=no

using Base.Iterators
using Logging
using Pipe: @pipe


"""
    compute_mpgl_sample_size

Parse the second line of the MPGL file, split out items by the space delimiter
" ", and count the number of entries. This row contains all the sample IDs with
genotype likelihoods represented in the remainder of the file.
"""
function compute_mpgl_sample_size(mpgl_path::String)

    # print a warning if the file doesn't end with .mpgl
    !endswith(mpgl_path, ".mpgl") && @warn "File does not have the mpgl suffix."

    # use a functional-programming-style pipeline to split out the sample size
    sample_size = @pipe readlines(mpgl_path) |>
        drop(_, 1) |>
        take(_, 1) |>
        collect |>
        split(_[1], " ") |>
        length

    # make sure the value is not zero
    @assert sample_size != 0 "Zero samples found. Aborting."

    return sample_size
end
precompile(compute_mpgl_sample_size, (String,))


"""
    write_starting_qs

Take the sample size and desired starting Q value, and write that Q-value on
one line per the number of samples. So, there are 60 samples, the file will have
60 lines containing that Q-value.
"""
function write_starting_qs(sample_size::Int64, starting_q::Float64, label::String)

    open("$(label)_starting_q.txt", "w") do writer
        _ = [println(writer, starting_q) for i in 1:sample_size];
    end

end
precompile(write_starting_qs, (Int64, Float64, String))


function main(args::Vector{String})

    # print warning if fewer than or more than 3 arguments were provided.
    length(args) != 3 && @warn """
    Three positional arguments must be provided:
    1. The relative or absolute path to an mpgl file for input to Entropy.
    2. A starting Q value to apply to all samples, e.g., 0.5.
    3. A label to be prepended onto the starting Q file name, e.g., "test"

    Proceeding with defaults for arguments 2 and 3.
    """

    # if 1 arg was provided, set the defaults
    if length(args) == 1
        (starting_q, label) = (0.5, "test")
    end

    # if more than one arg was provided, test it and set it
    if length(args) > 1
        @assert typeof(args[2]) == Float64 """
        The second argument must be a decimal between zero and one
        """
        starting_q = args[2]
    end

    # make sure the third arg follows the rules
    if length(args) > 2
        @assert typeof(args[3]) == String """
        The type of the third argument, the output file label, must be string.
        """
        label = args[3]
    end

    # run the functions
    @pipe compute_mpgl_sample_size(args[1]) |>
        write_starting_qs(_, starting_q, label)

end

# use short-circuiting to see if the function should run when the script is
# executed
length(ARGS) > 0 && isfile(ARGS[1]) && main(ARGS)
