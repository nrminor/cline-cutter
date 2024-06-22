#!/usr/bin/env -S julia --threads auto --history-file=no --startup-file=no --compiled-modules=no

using CSV, DataFrames, VariantCallFormat, VCFTools
using Base.Threads: @threads


"""
The struct `DatasetInfo` holds metadata about the VCF-formatted
dataset as a whole, including the number of markers (L), the sample
size (N), a 1-dimensional array of sample IDs, and the ploidy of
the species represented.
"""
struct DatasetInfo
    L::Int64
    N::Int64
    Sample_IDs::Vector{String}
    Ploidy::Int64
end


"""
The struct `EntropyHeader` is a simple container for the first and
and second header lines expected by Entropy in the MPGL format. The
first line is the sample size and the number of markers, and the second
line is the list of sample IDs. As with the rest of the file format,
all elements are space-delimited.
"""
struct EntropyHeader
    line1::String
    line2::String
end


"""
Take the VCF file and process it so that it can be rendered as
a Julia dataframe.
"""
function to_mpgl_df(vcf_file::String, vcf_meta::DatasetInfo)::DataFrame

    # create a sample data frame to fill with the correct PL values
    mpgl_df = DataFrame(fill("", (vcf_meta.L, vcf_meta.N)), :auto)
    rename!(mpgl_df, vcf_meta.Sample_IDs)

    # prepare to collect unique locus identifiers
    chrom = Vector{String}(undef, vcf_meta.L)
    mpgl_df.chrom = chrom

    # iterate through each record and extract the information needed for
    # mpgl format
    vcf_reader = VCF.Reader(open(vcf_file, "r"))
    for (i, record) in enumerate(vcf_reader)
        # fill in the current marker's identifier
        @inbounds mpgl_df.chrom[i] = VCF.chrom(record)

        # parse out genotypes, double checking that enough PL's are present
        genotypes = VCF.genotype(record, :, "PL")

        # fill in the correct genotype
        for (j, genotype) in enumerate(genotypes)
            @inbounds mpgl_df[i,j] = replace(genotype, "," => " ")
        end
    end

    # reorder so that chroms are first
    select!(mpgl_df, :chrom, :)

    return mpgl_df
end


"""
Collect information about the VCF so that a dataframe with MPGL
formatted information, along with the required header, can be
compiled.
"""
function vcf_to_mpgl(vcf_file::String, ploidy::Int64)

    # check that the VCF exists and log out
    @assert isfile(vcf_file) "Provided VCF path, '$vcf_file', does not exist in the working directory."
    println("Formatting input genotype likelihood data for the Entropy program from VCF $vcf_file")

    # parse out the sample IDs and the varying loci, along with counts
    sample_ids = sampleID(vcf_file)
    sample_size = length(sample_ids)
    n_loci = nrecords(vcf_file)

    # initialize the structs to bundle data
    vcf_meta = DatasetInfo(
        n_loci,
        sample_size,
        sample_ids,
        ploidy,
    )
    header = EntropyHeader(
        "$sample_size $n_loci",
        join(sample_ids, " ")
    )

    # get the Phred-scaled genotype likelihoods as a dataframe
    mpgl_df = to_mpgl_df(vcf_file, vcf_meta)

    return (header, mpgl_df)

end


"""
Write out the MPGL header and dataframe.
"""
function write_mpgl(header::EntropyHeader, mpgl_df::DataFrame, name_prefix::String)

    # construct the output name based on the input prefix
    tmp_name = "$name_prefix.tmp.mpgl"
    output_name = "$name_prefix.mpgl"

    # start with the header
    open(tmp_name, "w") do writer
        println(writer, header.line1)
        println(writer, header.line2)
    end

    # finish with the table data
    open(tmp_name, "a") do appender
        CSV.write(
            appender, mpgl_df;
            append=true, delim=' ', quotestrings=false, writeheader=false
        )
    end

    # remove quotes
    open(tmp_name, "r") do reader
        open(output_name, "w") do writer
            for line in eachline(reader)
                println(writer, replace(line, "\"" => ""))
            end
        end
    end
    rm(tmp_name)

end


"""
Script entrypoint
"""
function convert_all_files()

    # search for VCF files that can be converted
    workingdir = length(ARGS) == 1 ? ARGS[1] : "."
    vcf_files = [file for file in readdir(workingdir; join=true) if occursin(".vcf", file)]

    # double check that VCF files were actually found
    @assert length(vcf_files) > 0 "No VCF files found in the provided working directory"

    # sort by size so that the smallest comes first, thereby limiting first-
    # compilation time
    vcf_sizes = [filesize(file) for file in vcf_files]
    vcf_files = vcf_files[sortperm(vcf_sizes)]

    # process each available VCF
    @threads for vcf_file in vcf_files
        # convert to MPGL dataframe
        header, mpgl_df = vcf_to_mpgl(vcf_file, 2)

        # parse out the input file simple name
        simple_name = replace(vcf_file, ".vcf" => "")

        # write out the new file and remove
        write_mpgl(header, mpgl_df, simple_name)
    end

end

convert_all_files()
