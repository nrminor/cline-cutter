#!/usr/bin/env -S julia --threads auto --gcthreads 3

using CSV, DataFrames

const SAMPLING_REGIME::String = ARGS[1]
const FITTING_LOG::String = ARGS[2]

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
function parse_fitting_log(fitting_log::String, regime::String)::String

    @assert isfile(fitting_log) "Provided path to cline fitting R script log does not exist"

    open(fitting_log, "r") do reader
        # read the lines of the file
        lines = readlines(reader);

        # parse out the names of the models
        model_names = [replace(line,
            "Fitting model labeled '"=>"",
            "':"=>"",
        ) for line in lines if contains(line, "Fitting model labeled")];

        # make sure it was able to parse model names with the expected susbtrings
        @assert (
            length(model_names) > 0
        ) "Unable to parse any model names in the provided log $fitting_log"

        # parse our the Metropolis scores for each model
        model_scores = [replace(line,
            "The Metropolis acceptance rate was "=>"",
        ) for line in lines if contains(line, "The Metropolis acceptance rate was")]

        # make sure it was able to parse model scores with the expected susbtrings
        @assert (
            length(model_scores) > 0
        ) "Unable to parse any model Metropolis scores in the provided log $fitting_log"

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
            "Metropolis Acceptance Rate" => model_scores
        );
        sort!(score_df, ["Metropolis Acceptance Rate"], rev = true)
        CSV.write("$(regime)_scores.csv", score_df);
    end

    return "$(regime)_scores.csv"
end

const WRITEOUT_NAME::String = parse_fitting_log(
    FITTING_LOG,
    SAMPLING_REGIME
)
println("Parsed logging information written out to $WRITEOUT_NAME")
