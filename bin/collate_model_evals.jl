#!/usr/bin/env -S julia --threads auto --gcthreads=3 --history-file=no --startup-file=no

using CSV, DataFrames
using Pipe: @pipe

"""
Function `parse_fitting_log()` takes the peculiarly formatted logs produced by the
`hzar` R package in project script `cline_fitting.R` and extracts two pieces of
information from it:
1. The names provided for each competing cline model in `cline_fitting.R`.
2. The Metropolis acceptance rates for each model after run through an MCMC chain.

This information can be used to evaluate the performance of each model with respect
to their explorations of genomic cline parameter space, where higher Metropolis
algorithm acceptance rates are better.

**Note!** The approach used here is highly error prone because of a simplifying
assumption it makes: It assumes that *all model names are listed in the same order
as their corresponding model scores*.

Assertions are used to suss out possible mismatches in the numbers of names and scores,
but the user is expected let `cline_fitting.R` work as written so that names always
correspond with scores. This approach may be replaced in the future with iterative
construction of structs that will make it clearer when a name or score is missing, e.g.:
```
struct ModelScores
    name::String,
    score::Float32
end
````
"""
function parse_fitting_log(fitting_log::String, regime::String)::Tuple{DataFrame,String}

    @assert isfile(fitting_log) "Provided path to cline fitting R script log does not exist"

    output_name = "$(regime)_regime_metropolis_rates.csv"

    open(fitting_log, "r") do reader
        # read the lines of the file
        lines = readlines(reader)

        # parse out the names of the models
        special_char_pattern = r"[!-,\.-\/:-@\[-\^`\{-~\s]+"
        model_names = [
            replace(
                line,
                "Fitting model labeled '" => "",
                "':" => "",
                special_char_pattern => "",
            ) for line in lines if contains(line, "Fitting model labeled")
        ]

        # make sure it was able to parse model names with the expected susbtrings
        @assert (length(model_names) > 0) "Unable to parse any model names in the provided log $fitting_log"

        # parse our the Metropolis scores for each model
        model_scores = [
            replace(line, "The Metropolis acceptance rate was " => "") for
            line in lines if contains(line, "The Metropolis acceptance rate was")
        ]

        # make sure it was able to parse model scores with the expected susbtrings
        @assert (length(model_scores) > 0) "Unable to parse any model Metropolis scores in the provided log $fitting_log"

        # make sure that there are the same numbers of scores and names
        @assert (
            length(model_names) == length(model_scores)
            ) """
            Different numbers of model names and model acceptance rates detected:
            Model names: $model_names
            Model scores: $model_scores
            Please double check that the provided log contains a score and a name
            for each model.
            """

        # structure as a dataframe and write it out
        score_df = DataFrame(
            "Model Label" => model_names,
            "Metropolis Acceptance Rate" => model_scores,
        )
        sort!(score_df, ["Metropolis Acceptance Rate"], rev = true)
        CSV.write(output_name, score_df)

        return (score_df, output_name)
    end
end;
precompile(parse_fitting_log, (String, String))

"""
The function `join_model_evals()` uses a simple leftjoin to combine AIC information
and Metropolis acceptance rates into one table listing both statistics for each
model fitted in the R script `cline_fitting.R`.
"""
function join_model_evals(score_df::DataFrame, aic_file::String, sampling_regime::String)

    aic_df = CSV.read(aic_file, DataFrame)
    writeout_name = "$(sampling_regime)_regime_model_evals.csv"

    @assert "model" in names(aic_df) "Expected 'model' column is missing in the AIC score table."
    @assert "AICc" in names(aic_df) "Expected 'AICc' column is missing in the AIC score table."

    joined_evals =
        @pipe score_df |> leftjoin(_, aic_df, on = Symbol("Model Label") => :model)

    CSV.write(writeout_name, joined_evals)

    return writeout_name

end;
precompile(join_model_evals, (DataFrame, String, String))

"""
Script entrypoint
"""
function main()

    sampling_regime = ARGS[1]
    fitting_log = ARGS[2]

    score_df, writeout_name = parse_fitting_log(fitting_log, sampling_regime)
    println("Parsed logging information written out to $writeout_name")

    aic_files = [
        file for
        file in readdir() if contains(file, "aic.tsv") & contains(file, sampling_regime)
    ]
    if length(aic_files) == 1
        full_table = join_model_evals(score_df, aic_files[1], sampling_regime)
        println(
            """
            Model Akaike Information Criteria and Metropolis acceptance rates written
            to $full_table in the current working directory. Deleting input tables.
            """
            )
        rm(writeout_name)
        rm(aic_files[1])
    end

end;
precompile(main, ())

length(ARGS) > 0 && main()
