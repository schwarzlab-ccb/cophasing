"""
Build the genome-wide bin table that `cooler load` matches contacts against.

Extracted from 04_generate_cool_file_permutation_test_results.py, where every
chromosome task used to build this file itself and serialise on a lock file.
Running it once as its own step removes both the lock and the repeated download.

Code by Claudia Robens
"""

import argparse

import pandas as pd
import requests

UCSC_CHROM_SIZES = "https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"

parser = argparse.ArgumentParser(description='Write a BED of fixed-size genomic windows')
parser.add_argument('--resolution', type=int, required=True,
                    help='Bin size in bp, e.g. 40000')
parser.add_argument('--out', type=str, required=True,
                    help='Path of the BED file to write')
args = parser.parse_args()

response = requests.get(UCSC_CHROM_SIZES, timeout=30)
response.raise_for_status()

hg38_chrom_sizes = {}
for line in response.text.strip().split("\n"):
    chrom, size = line.split("\t")
    hg38_chrom_sizes[chrom] = int(size)

autosomes = [f"chr{i}" for i in range(1, 23)]
chromosome_sizes = {c: hg38_chrom_sizes[c] for c in autosomes}

bed_data = []
for chrom, size in chromosome_sizes.items():
    start = 0
    while start < size:
        # The final window of a chromosome is truncated at the chromosome end
        end = min(start + args.resolution, size)
        bed_data.append([chrom, start, end])
        start += args.resolution

bed_df = pd.DataFrame(bed_data, columns=['Chromosome', 'start', 'end'])
bed_df.to_csv(args.out, sep='\t', index=False, header=False)
print(f"{len(bed_df)} genomic windows at resolution {args.resolution} written to {args.out}")
