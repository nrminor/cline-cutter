#!/usr/bin/env python3

"""
DOWNSAMPLE A MULTI-SAMPLE VCF BY GEOGRAPHIC LOCATION
----------------------------------------------------

This script takes a multi-sample VCF, cross-references its sample names
and latitude/longitude coordinates with a metadata file, and downsamples
each VCF according to its population "cluster". To do so, it computes a
geographic distance matrix for all samples according to their coordinates,
and groups samples together by their proximity. It then downsamples the
VCFs randomly, evenly, and unevenly across population clusters (more 
granular explanation in the function docstrings below). The purpose of
this downsampling is to stress-test population genomic analyses whose
resolution or precision may be affected not only by sample size, but
also the density of samples across space.
"""

import os
import argparse
import subprocess
from subprocess import Popen
import glob
from math import ceil
import gzip
# import random
import polars as pl
import numpy as np
from geopy.distance import distance
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

def parse_command_line_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", "-m",
                        type=str,
                        required=True,
                        help="Metadata file with sample IDs and decimal latitudes and longitudes.")
    parser.add_argument("--distance_threshold", "-d",
                        type=int,
                        default=100,
                        required=False,
                        help="The maximum kilometer distance at which a new cluster \
                            should be assigned.")
    parser.add_argument("--proportion", "-p",
                        type=float,
                        default=0.5,
                        required=False,
                        help="Proportion of samples to retain in downsampling.")
    parser.add_argument("--vcf", "-v",
                        type=str,
                        required=True,
                        help="Multisample VCF to subset.")
    parser.add_argument("--cores", "-c",
                        type=int,
                        default=3,
                        required=False,
                        help="Number of cores to use in multiprocessing.")
    parser.add_argument("--seed", "-s",
                        type=int,
                        default=14,
                        required=False,
                        help="Seed for random number generation.")
    args = parser.parse_args()
    return args.metadata, args.distance_threshold, args.proportion, args.seed, args.vcf

def unify_sample_names(metadata_df: pl.DataFrame, vcf_path: str) -> str:
    """
    This function double checks the VCF at the provided path to make
    sure the sample names in the VCF are the same as the sample names
    in the metadata. If they aren't, it looks for the metadata sample
    names as substrings of the VCF sample names. If all metadata 
    sample names are in fact substrings of the VCF sample names, it
    creates a tab-delimited text file to guide BCFTools, which it
    uses in a subprocess to rename VCF sample names.

    Args:
    metadata_df - A Polars DataFrame of the sample metadata.
    vcf_path - a character string filepath leading to the input VCF.

    Returns:
    A character string path to the new VCF output by VCFtools.
    """

    # create a list of sample names from the VCF
    vcf_samples = []
    if vcf_path.endswith(".gz"):
        with gzip.open(vcf_path, 'rb') as vcf_file:
            for line in vcf_file:
                if line.startswith(b'#CHROM'):
                    # Split the line by tabs, and extract the sample IDs
                    # starting from the 9th column (0-indexed)
                    vcf_samples = line.decode('utf-8').strip().split("\t")[9:]
                    break
    else:
        with open(vcf_path, 'r', encoding="utf-8") as vcf_file:
            for line in vcf_file:
                if line.startswith('#CHROM'):
                    # Split the line by tabs, and extract the sample IDs
                    # starting from the 9th column (0-indexed)
                    vcf_samples = line.strip().split('\t')[9:]
                    break

    # create a list of the sample names from the metadata
    meta_samples = metadata_df.select(
        "Sample ID"
    ).drop_nulls().to_series().to_list()

    # match the two lists, permitting no samples to be unmatched
    matched_sample_ids = []
    for meta_sample in meta_samples:
        for vcf_sample in vcf_samples:
            if meta_sample in vcf_sample:
                matched_sample_ids.append((vcf_sample, meta_sample))
                break

    # write a two column text file, where the first column is names
    # in the VCF (old names), and the second column is names in the
    # metadata (new names)
    with open('new_sample_names.txt', 'w', encoding="utf-8") as out_handle:
        for line in matched_sample_ids:
            print(f"{line[0]}\t{line[1]}", file=out_handle)

    # create the BCFTools command and run it in a subprocess
    new_vcf = "renamed_samples.vcf.gz"
    with open(new_vcf, "w", encoding="utf-8") as vcf_handle:
        process = subprocess.Popen(("bcftools", "reheader", "--samples",
                           "new_sample_names.txt", vcf_path), 
                           stdout=subprocess.PIPE)
        subprocess.run(["gzip", "-c"], check=True, stdin=process.stdout, stdout=vcf_handle)
        process.wait()

    # clear out temporary sample name file
    os.remove("new_sample_names.txt")

    return new_vcf

def assign_clusters(metadata_df: pl.DataFrame, max_distance: int) -> pl.DataFrame:
    """
    This function uses a Vincenty distance matrix to assign a
    cluster index to each sample, which will be used for downsampling
    downstream.

    Args:
    metadata_df - a Polars DataFrame that includes a column of Sample
    IDs, a decimal latitude column, and a decimal longitude column.

    Returns:
    A Polars DataFrame with a Sample ID column and a cluster index
    column, listing which samples land in which clusters.
    """

    # select down to sample ID, longitude, and latitude, sort
    # alphabetically by sample ID
    metadata_df = metadata_df.select(
        pl.col(["Sample ID", "Latitude", "Longitude"])
    ).sort(by="Sample ID")

    # create a dictionary of zipped latitude, longitude values called 'coords'
    sample_ids = metadata_df.select(pl.col("Sample ID")).to_series().to_list()
    lats = metadata_df.select(pl.col("Latitude")).to_series().to_list()
    longs = metadata_df.select(pl.col("Longitude")).to_series().to_list()
    coordinates = list(zip(lats, longs))

    # create an empty array to iterate through
    distmat = np.zeros(shape=(len(sample_ids), len(sample_ids)))

    # iterate through each sample ID, compute vincenty distances to all
    # sample IDs, and replace a "column" of the array with the resulting list
    for i, coordinate in enumerate(coordinates):
        distmat[::,i] = [distance(coordinate,
                                  other_coord,
                                  ellipsoid='WGS-84').kilometers for other_coord in coordinates]

    # cluster the distances using hierarchical clustering with complete linkage
    linked = linkage(squareform(distmat), method='complete')

    # Form flat clusters based on the threshold
    clusters = fcluster(linked, max_distance, criterion='distance')

    # Assign cluster indices based on distances
    cluster_assignments = pl.DataFrame({
        "sample": sample_ids,
        "cluster": clusters
    })

    return cluster_assignments

def create_sample_lists(metadata_df: pl.DataFrame,
                        max_distance: int,
                        proportion: float,
                        seed: int) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """
    This function takes the outputs from function assign_clusters
    and uses it to downsample the list of Sample IDs by three
    regimes: 1) random sampling that is unsensitive to sample 
    location; 2) even sampling that downsamples all population
    clusters evenly; And 3) uneven sampling, which samples some
    clusters more than others. Note that the function deals with
    singleton clusters by simply keeping them and not downsampling.
    This may be replaced with a `random.choice()` call to determine
    whether a singleton is kept or not in the future.

    Args:
    metadata_df - A Polars DataFrame of Sample IDs and coordinates.

    Returns:
    Three polars dataframes, each of which lists sample IDs for
    one of the above downsampling regimes. These lists can be
    used to downsample a multisample VCF with the same sample IDs.
    """

    # retrieve cluster dataframe
    cluster_df = assign_clusters(metadata_df, max_distance)

    # random sampling
    random_df = cluster_df.sample(fraction=proportion, seed=seed)
    random_df = random_df.select(pl.col("sample"))

    # even sampling
    even_df = cluster_df.clear()
    even_dict = cluster_df.partition_by(by="cluster", as_dict=True)
    for key in even_dict:
        sub_df = even_dict.get(key)
        if sub_df.shape[0] == 1:
            # whether_to_keep = bool(random.getrandbits(1))
            even_df.vstack(sub_df, in_place=True)
        else:
            sub_df = sub_df.sample(fraction=proportion, seed=seed)
            even_df.vstack(sub_df, in_place=True)
    even_df = even_df.select(pl.col("sample"))

    # uneven sampling
    clusters = cluster_df.select("cluster").unique(subset="cluster").to_series().to_list()
    uneven_down = np.random.choice(clusters, size=ceil(len(clusters) / 3))
    unchanged = cluster_df.filter(
        pl.col("cluster") != uneven_down[0]
    )
    uneven_df = cluster_df.filter(
        pl.col("cluster") == uneven_down[0]
    ).sample(
        fraction=proportion, seed=seed
    ).vstack(unchanged).select(pl.col("sample"))

    return (random_df, even_df, uneven_df)

def main():
    """
    Main handles file I/O and a VCFtools subprocess to bring
    the above functions together.

    Args:
    metadata_path - Path to the metadata spreadsheet with sample IDs, decimal 
    latitude, and decimal longitude.
    max_distance - The maximum number of kilometers for two samples to be 
    placed in the same cluster.
    vcf_path - Path to the multisample VCF to downsample based on geographic
    clusters.

    Returns:
    None
    """

    # parse arguments from the command line
    metadata_path, max_distance, proportion, seed, vcf_path = parse_command_line_args()

    # read metadata into a Polars DataFrame
    metadata = pl.read_excel(metadata_path,
                                read_csv_options={"null_values": ["NA", ""]}).select(
        pl.col(["Sample ID", "Latitude", "Longitude"])
    ).drop_nulls()

    # create sample lists with which to filter the VCF
    random_df, even_df, uneven_df = create_sample_lists(metadata, max_distance, proportion, seed)

    # make sure sample names in the VCF are the same as the metadata
    new_vcf = unify_sample_names(metadata, vcf_path)

    # write out sample lists
    sample_types = ["random_sample", "even_sample", "uneven_sample"]
    random_df.write_csv(f"{sample_types[0]}.txt", separator="\t", has_header=False)
    even_df.write_csv(f"{sample_types[1]}.txt", separator="\t", has_header=False)
    uneven_df.write_csv(f"{sample_types[2]}.txt", separator="\t", has_header=False)

    # prepare commands for parallel filtering
    cmds_list = [["vcftools", "--gzvcf", new_vcf, "--keep",
                  f"{sample_type}.txt", "--recode", "--out", 
                  sample_type] for sample_type in sample_types]

    # use Popen to run filtering in parallel
    processs_list = [Popen(cmd) for cmd in cmds_list]
    for process in processs_list:
        process.wait()

    # remove VCFtools log files and other intermediate files
    os.remove("renamed_samples.vcf.gz")
    for file in glob.glob("*.log"):
        os.remove(file)


if __name__ == "__main__":
    main()
