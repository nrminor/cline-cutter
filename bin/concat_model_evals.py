#!/usr/bin/env python3

"""
Usage:
```
python3 concat_model_evals.py <PATTERN>
```
"""
import os
import sys
from typing import List, Optional, Tuple

import polars as pl


def find_files(match_pattern: str) -> Tuple[str]:
    """
    Collect all files that match a certain pattern and return in a tuple.
    """
    cwd_files = os.listdir(".")
    eval_files = tuple(
        file for file in cwd_files if match_pattern in file and os.path.isfile(file)
    )

    return eval_files


def check_files(
    eval_files: Tuple[str],
    checked_files: Optional[List[str]] = None,
    remaining_files: Optional[int] = None,
) -> List[str]:
    """
    Recursively check files for consistency in column numbers, names, and types.
    """

    # initialize checked files if this is the first call of the function
    checked_files = [] if checked_files is None else checked_files

    # if on the first encountered file, which is to say the file at the end
    # of the list of files to check, make sure the file can be read, but otherwise
    # effectively skip and begin recursion
    if remaining_files is None:
        _ = pl.read_csv(eval_files[-1], n_rows=1)
        return check_files(
            eval_files, checked_files.append(eval_files[-1]), len(eval_files) - 1
        )

    # if on any file other than the one at the end of the list of files to check,
    # set the index for the previous file we'll use for comparison
    prev_idx = remaining_files + 1

    # read the files to compare, taking one row from each
    file_to_check = pl.read_csv(eval_files[remaining_files], n_rows=1)
    previous_passing = pl.read_csv(eval_files[prev_idx], n_rows=1)

    # Run assertions to check for column numbers, names, and data types
    assert (
        previous_passing.shape[1] == file_to_check.shape[1]
    ), f"Column number for {eval_files[remaining_files]} does not match {eval_files[prev_idx]}"
    assert all(
        col in previous_passing.columns for col in file_to_check.columns
    ), f"Column names in {eval_files[remaining_files]} do not match {eval_files[prev_idx]}"
    assert (
        file_to_check.schema == previous_passing.schema
    ), f"Data types in {eval_files[remaining_files]} do not match {eval_files[prev_idx]}"

    # if all is well and the current file passed, append it to the checked files list
    checked_files.append(eval_files[remaining_files])

    # early return if on the final file. Otherwise, recurse.
    if remaining_files == 0:
        return checked_files
    return check_files(eval_files, checked_files, remaining_files - 1)


def main() -> None:
    """
    Script entrypoint
    """

    # pull in the first argument as the pattern to use for finding files
    match_pattern = sys.argv[1]

    # find the files of model evaluations
    eval_files = find_files(match_pattern)

    # make sure files were actually found
    assert (
        len(eval_files) > 0
    ), "No model evaluation files found in the current working directory."

    # check that all the files have the same columns and types
    checked_files = check_files(eval_files)

    # read the files into an immutable tuple
    eval_df_list = (pl.scan_csv(file) for file in checked_files)

    # append each file into a single dataframe and write out
    pl.concat(eval_df_list).sink_csv("concatenated_evals.csv")


if __name__ == "__main__":
    main()
