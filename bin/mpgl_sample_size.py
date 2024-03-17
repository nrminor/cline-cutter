#!/usr/bin/env python3

"""
```
usage: mpgl_sample_size.py [-h] --mpgl MPGL [--starting_q STARTING_Q] [--label LABEL]

options:
  -h, --help            show this help message and exit
  --mpgl MPGL, -m MPGL  The relative or absolute path to an mpgl file for input to Entropy.
  --starting_q STARTING_Q, -q STARTING_Q
                        A starting Q value to apply to all samples, e.g., 0.5.
  --label LABEL, -l LABEL
                        A label to be prepended onto the starting Q file name, e.g., 'test'
```
"""

import argparse
import asyncio
from os.path import isfile
from pathlib import Path
from typing import Tuple

import polars as pl


def parse_command_line_args() -> Tuple[Path, float, str]:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mpgl",
        "-m",
        type=Path,
        required=True,
        help="The relative or absolute path to an mpgl file for input to Entropy.",
    )
    parser.add_argument(
        "--starting_q",
        "-q",
        type=float,
        required=False,
        default=0.5,
        help="A starting Q value to apply to all samples, e.g., 0.5.",
    )
    parser.add_argument(
        "--label",
        "-l",
        type=str,
        required=False,
        default="test",
        help="A label to be prepended onto the starting Q file name, e.g., 'test'",
    )
    args = parser.parse_args()

    return (args.mpgl, args.starting_q, args.label)


async def compute_mpgl_sample_size(mpgl_path: Path):
    """
        compute_mpgl_sample_size

    Parse the second line of the MPGL file, split out items by the space delimiter
    " ", and count the number of entries. This row contains all the sample IDs with
    genotype likelihoods represented in the remainder of the file.
    """

    # print a warning if the file doesn't end with .mpgl
    if not str(mpgl_path).endswith(".mpgl"):
        print(
            "File does not have the mpgl suffix. Successful parsing can't be guaranteed."
        )

    # use a functional-programming-style pipeline to split out the sample size
    sample_size: int
    with open(mpgl_path, "r", encoding="utf8") as file:
        sample_size = len(
            pl.Series("lines", list(file.readlines()))
            .gather(1)
            .str.split(" ")
            .to_list()[0]
        )

    assert sample_size != 0, "Zero samples found in input file. Aborting."

    return sample_size


async def write_starting_qs(sample_size: int, starting_q: float, label: str) -> str:
    """
        write_starting_qs

    Take the sample size and desired starting Q value, and write that Q-value on
    one line per the number of samples. So, there are 60 samples, the file will have
    60 lines containing that Q-value.
    """

    out_name = f"{label}_starting_q.txt"

    with open(out_name, "w", encoding="utf8") as out_handle:
        out_handle.writelines([f"{str(starting_q)}\n" for _ in range(sample_size)])

    return out_name


async def main() -> None:
    """
    Main runs the above functions in an asynchronous runtime.
    """
    mpgl, starting_q, label = parse_command_line_args()

    assert isfile(
        mpgl
    ), "Provided MPGL file path does not point toward a file that exists on the current filesystem."

    sample_size = await compute_mpgl_sample_size(mpgl)

    out_name = await write_starting_qs(sample_size, starting_q, label)

    print(f"Starting Q's written to file {out_name}")


if __name__ == "__main__":
    asyncio.run(main())
