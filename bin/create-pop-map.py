#!/usr/bin/env python3

"""
CREATE TAB-DELIMITED POPULATION MAP FILE
----------------------------------------
"""

import argparse
import polars

def parse_command_line_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--samplesheet", "-s",
                        type=str,
                        required=True,
                        help="Metadata file with sample IDs and taxon shorthands or \
                            population codes. Note that these sample IDs must match \
                            the sample IDs in the VCFs or whatever other population \
                            genomic format you are using.")
    parser.add_argument("--project_name", "-p",
                       type=str,
                       required=False,
                       default="population_map",
                       help="what to name the population map fie.")
    args = parser.parse_args()
    return args.samplesheet, args.project_name

def main():
    """
    Main handles file I/O and samplesheet subsetting to create a population map.
    """

    # parse arguments from the command line
    samplesheet, project_name = parse_command_line_args()

    # read samplesheet into memory
    metadata = polars.read_excel(samplesheet,
                                read_csv_options={"null_values": ["NA", ""]})

    # make sure the required column names are in that samplesheet
    assert "Sample ID" in metadata.columns, "A column named 'Sample ID' must \
        be present in the samplesheet."
    assert "Taxon shorthand" in metadata.columns, "A column named 'Taxon shorthand' \
        must be present in the samplesheet."

    # filter to the columns we care about
    pop_map = metadata.select(
        polars.col(["Sample ID", "Taxon shorthand"])
    ).drop_nulls()

    # write new pop map with ".pop" appendage
    pop_map.write_csv(f"{project_name}.pop", separator="\t", has_header=False)


if __name__ == "__main__":
    main()
