#!/usr/bin/env python3

"""
usage:
```
python3 collate_model_evals.py <SUBSAMPLING_LABEL> <MODELING_LOG>
```
"""

import os
import re
import sys
from pathlib import Path
from typing import Optional, Tuple
from dataclasses import dataclass

import polars as pl


@dataclass(frozen=True)
class ReplicateData:
    """
    ReplicateData stores metadata about the cline-cutter run, namely
    the downsampling regime used, the proportion of samples downsampled,
    and the random number seed used for this replicate.
    """

    regime: str
    proportion: Optional[float]
    seed: Optional[int]


def parse_fitting_log(fitting_log: Path, regime: str) -> Tuple[pl.LazyFrame, str]:
    """
    Function `parse_fitting_log()` takes the peculiarly formatted logs produced by the
    `hzar` R package in project script `cline_fitting.R` and extracts two pieces of
    information from it:
    1. The names provided for each competing cline model in `cline_fitting.R`.
    2. The cline center and width estimates from each model.
    3. The Metropolis acceptance rates for each model after run through an MCMC chain.

    This information can be used to evaluate the performance of each model with respect
    to their explorations of genomic cline parameter space, where higher Metropolis
    algorithm acceptance rates are better.

    **Note!** The approach used here is highly error prone because of a simplifying
    assumption it makes: It assumes that *all model names are listed in the same order
    as their corresponding model scores*.

    Assertions are used to suss out possible mismatches in the numbers of names and scores,
    but the user is expected let `cline_fitting.R` work as written so that names always
    correspond with scores. This approach may be replaced in the future with iterative
    construction of dataclasses that will make it clearer when a name or score is missing, e.g.:
    ```
    @dataclass
    class ModelScores:
        name: str,
        score: float,
    ````
    """

    assert os.path.isfile(
        fitting_log
    ), "Provided path to cline fitting R script log does not exist"

    output_name = f"{regime}_regime_metropolis_rates.csv"

    # collect the lines of the log file
    with open(fitting_log, "r", encoding="utf8") as log_handle:
        lines = log_handle.readlines()

    # parse out the names of the models
    special_char_pattern = r"[!-,\.-\/:-@\[-\^`\{-~\n\r]+"
    cleaned_lines = [re.sub(special_char_pattern, "", line) for line in lines]
    model_names = [
        (line.replace("Fitting model labeled ", "").replace("':", ""))
        for line in cleaned_lines
        if "Fitting model labeled" in line
    ]

    # make sure it was able to parse model names with the expected susbtrings
    assert (
        len(model_names) > 0
    ), f"Unable to parse any model names in the provided log {fitting_log}"

    # parse out the cline centers and widths
    cline_center = [
        (
            line.replace(
                f'[1] "Estimated cline center for the {regime} downsampling regime: ',
                "",
            ).replace("':", "")
        )
        for line in cleaned_lines
        if "Estimated cline center for the" in line
    ]
    cline_width = [
        (
            line.replace(
                f'[1] "Estimated cline width for the {regime} downsampling regime: ', ""
            ).replace("':", "")
        )
        for line in cleaned_lines
        if "Estimated cline width for the" in line
    ]

    # make sure it was able to parse cline statistics with the expected susbtrings
    assert (
        len(cline_center) > 0
    ), f"Unable to parse any cline centers in the provided log {fitting_log}"
    assert (
        len(cline_width) > 0
    ), f"Unable to parse any cline widths in the provided log {fitting_log}"

    # parse our the Metropolis scores for each model
    model_scores = [
        float(line.replace("The Metropolis acceptance rate was ", ""))
        for line in lines
        if "The Metropolis acceptance rate was" in line
    ]

    # make sure it was able to parse model scores with the expected susbtrings
    assert (
        len(model_scores) > 0
    ), f"Unable to parse any model Metropolis scores in the provided log {fitting_log}"

    # make sure that there are the same numbers of scores and names
    assert len(model_names) == len(model_scores), f"""
        Different numbers of model names and model acceptance rates detected:
        Model names: {model_names}
        Model scores: {model_scores}
        Please double check that the provided log contains a score and a name
        for each model.
        """

    # structure as a dataframe
    score_df = pl.LazyFrame(
        {"Model Label": model_names, "Metropolis Acceptance Rate": model_scores}
    )

    # pull in cline statistics with a join
    cline_df = score_df.filter(
        pl.col("Model Label").str.contains("free")
        & pl.col("Model Label").str.contains("none")
    ).with_columns(
        pl.lit(cline_width).alias("Cline Width Estimate"),
        pl.lit(cline_center).alias("Cline Center Estimate"),
    )
    (
        score_df.join(cline_df, on="Model Label", how="left")
        .sort("Metropolis Acceptance Rate", descending=True)
        .sink_csv(output_name)
    )

    return (score_df, output_name)


def join_model_evals(score_df: pl.LazyFrame, aic_file: Path, rep_data: ReplicateData):
    """
    The function `join_model_evals()` uses a simple leftjoin to combine AIC information
    and Metropolis acceptance rates into one table listing both statistics for each
    model fitted in the R script `cline_fitting.R`.
    """

    aic_df = pl.scan_csv(aic_file, separator="\t")
    writeout_name = f"{rep_data.regime}_regime_p{rep_data.proportion}_s{rep_data.seed}_model_evals.csv"

    assert (
        "model" in aic_df.columns
    ), "Expected 'model' column is missing in the AIC score table."
    assert (
        "AICc" in aic_df.columns
    ), "Expected 'AICc' column is missing in the AIC score table."

    (
        aic_df.join(score_df, how="left", left_on="model", right_on="Model Label")
        .with_columns(
            pl.col("model")
            .str.replace("cline_", "")
            .str.replace("Model", "")
            .str.split("_")
            .alias("model")
        )
        .with_columns(
            [
                pl.col("model").list.first().alias("Model bounding"),
                pl.col("model").list.last().alias("Cline side"),
                pl.lit(rep_data.regime).alias("Regime"),
                pl.lit(rep_data.proportion).alias("Proportion downsampled"),
                pl.lit(rep_data.seed).alias("Replicate seed"),
            ]
        )
        .drop("model")
        .select(
            [
                "Model bounding",
                "Cline side",
                "Regime",
                "Proportion downsampled",
                "Replicate seed",
                "Cline Center Estimate",
                "Cline Width Estimate",
                "AICc",
                "Metropolis Acceptance Rate",
            ]
        )
        .sort(["AICc", "Metropolis Acceptance Rate"], descending=True)
        .sink_csv(writeout_name)
    )

    return writeout_name


def main() -> None:
    """
    Script entrypoint
    """
    sampling_regime = sys.argv[1]
    fitting_log = sys.argv[2]

    # split out sampling regime, proportion, and seed and package as a portable dataclass
    split_id = sampling_regime.split("_")
    if len(split_id) == 3:
        rep_data = ReplicateData(
            split_id[0],
            split_id[1],
            split_id[2],
        )
    else:
        rep_data = ReplicateData(split_id[0], None, None)

    score_df, writeout_name = parse_fitting_log(fitting_log, sampling_regime)
    print(f"Parsed logging information written out to {writeout_name}")

    aic_files = [
        str(os.path.realpath(file))
        for file in os.listdir(".")
        if "aic.tsv" in str(os.path.realpath(file))
        and sampling_regime.split("_")[0] in str(os.path.realpath(file))
    ]

    if not len(aic_files) == 1:
        print("AIC file not found. Finishing here.")
        sys.exit(0)

    full_table = join_model_evals(score_df, aic_files[0], rep_data)
    assert os.path.isfile(full_table), f"Writing of {full_table} failed."
    print(
        f"""
        Model Akaike Information Criteria and Metropolis acceptance rates written
        to {full_table} in the current working directory. Deleting input tables.
        """
    )
    os.remove(writeout_name)


if __name__ == "__main__":
    main()
