"""
SCRIPT TO GENERATE COOL FILE OF NPMI MATRICES TO PLOT IN PYGENOMETRACKS
Script to keep NaN values and in a second option, replace them with 5s
Save both cool files in the same directory
Code by Claudia Robens
"""

import logging
logging.basicConfig(level=logging.INFO)  # Set the logging level as needed
logging.info("Script started.")


"""
Imports
"""
from post_permutation_test import find_and_load_pkl, process_permutation_results, process_permutation_results_thresh
from generate_npmi import coordinates_from_location_string, get_region_from_location_string, calculate_NPMI, get_ticks, get_npmi_scale, read_segregation_table


import pandas as pd
import numpy as np
from pathlib import Path
import gzip
import subprocess
import shlex
import fcntl
import os
import glob
import argparse
import requests


logging.info("imports done.")
print("\n")


# define args
"""
Command-line arguments
"""
parser = argparse.ArgumentParser()
parser.add_argument('--chr', type=str, required=True, help='chr number')
parser.add_argument('--output_dir', type=str, required=True,
                    help='Base output directory containing all pipeline subdirectories')
parser.add_argument('--cutoff', type=str, required=True,
                    help='The cutoff value which was used for Co-Phasing pipeline, e.g. "10Mb"')
parser.add_argument('--resolution', type=int, required=True,
                    help='The resolution which was used for Co-Phasing pipeline, e.g. 40000')
parser.add_argument('--gaussian_kernel_size', type=int, required=True,
                    help='Kernel size used for Gaussian filter')

args = parser.parse_args()
logging.info("args done.")

curated_dir = str(Path(args.output_dir) / '01_curated_segregation_tables')
permutation_test_dir = str(Path(args.output_dir) / '02_permutation_test')
bins_rm_file = str(Path(args.output_dir) / '01_curated_segregation_tables' / 'bins_rm_hap1_hap2.tsv')
thresh_ident = str(Path(args.output_dir) / '03_thresholds' / 'knee_point_table_threshold_steps.tsv')


# Print arguments to stdout
print("\n" + "="*40)
print("Arguments applied:")
for arg, value in vars(args).items():
    print(f"{arg}: {value}")
print("="*40 + "\n")


segregation_table_hap1 = read_segregation_table(Hap="hap1", file_path=curated_dir,
                                                 resolution=args.resolution, cutoff=args.cutoff)
segregation_table_hap2 = read_segregation_table(Hap="hap2", file_path=curated_dir,
                                                 resolution=args.resolution, cutoff=args.cutoff)
segregation_table_both = read_segregation_table(Hap="both", file_path=curated_dir,
                                                 resolution=args.resolution, cutoff=args.cutoff)

"""
calcualte NPMI matrices
"""
# Hap1
npmi_matrix_hap1 = calculate_NPMI (segregation_table_hap1, args.chr)
# Hap2
npmi_matrix_hap2 = calculate_NPMI (segregation_table_hap2, args.chr)
# Both
npmi_matrix_both = calculate_NPMI (segregation_table_both, args.chr)

subset_segtable_hap1 = get_region_from_location_string (segregation_table_hap1, args.chr)

print("="*40 + "\n")


"""
Read in permutation results
"""


perm_results = find_and_load_pkl(directory_path=permutation_test_dir, chrom=args.chr)
filtered_array = process_permutation_results_thresh(perm_results=perm_results, chrom=args.chr,
                                                     thresh_table=thresh_ident, bins_rm_f=bins_rm_file,
                                                     subset_segtable_hap1=subset_segtable_hap1,
                                                     gaussian_kernel_size=args.gaussian_kernel_size)


p_value_npmi_2sided_updated = filtered_array['p_value_npmi_2sided_updated']
p_value_log_filter_npmi = filtered_array['p_value_log_filter_npmi']
p_value_log_filter_thresh_npmi = filtered_array['p_value_log_filter_thresh_npmi']
p_value_log_filter_thresh_npmi_directionality = filtered_array['p_value_log_filter_thresh_npmi_directionality']


"""remove bins in NPMI (unphased and phased) based on WDF curation"""
bins_rm = pd.read_csv(bins_rm_file, sep='\t', index_col=0)
bins_rm_filtered = bins_rm[bins_rm['chrom'] == args.chr]
tbr = subset_segtable_hap1.index.get_level_values('start').isin(bins_rm_filtered['start'])  # to be removed


# Set columns and rows to NaN which are identified in WDF curation step
npmi_matrix_hap1[:, np.array(tbr)] = np.nan  # set columns to NaN
npmi_matrix_hap1[np.array(tbr), :] = np.nan  # set rows to NaN

npmi_matrix_hap2[:, np.array(tbr)] = np.nan  # set columns to NaN
npmi_matrix_hap2[np.array(tbr), :] = np.nan  # set rows to NaN

npmi_matrix_both[:, np.array(tbr)] = np.nan  # set columns to NaN
npmi_matrix_both[np.array(tbr), :] = np.nan  # set rows to NaN

"""
Generate dataframe with 'chr:start-end' information as rownames and column names
"""
def add_chr_column(array, chr, resolution):
    # Create an array with the specified strings
    chr_column = np.array([f'{chr}:{i}-{i+int(resolution)}' for i in range(0, len(array)*int(resolution), int(resolution))])

    # Add the new column to the beginning of the existing array
    array_with_chr_column = np.column_stack((chr_column, array))
    
    # Generate column and row names
    column_names = np.array([f'{chr}:{i}-{i+int(resolution)}' for i in range(0, array.shape[1]*int(resolution), int(resolution))])
    
    # Convert to Pandas DataFrame and set row names as index
    df = pd.DataFrame(data=array, columns=column_names)
    df.index = chr_column   
    
    return df


def create_long_matrix(df):
    # Resetting the index to include row names as a regular column
    dataframe_reset = df.reset_index()

    # Melt the DataFrame to long format
    long_matrix = pd.melt(dataframe_reset, id_vars=['index'], var_name='Column', value_name='value')

    # Rename the columns
    long_matrix.columns = ['Row', 'Column', 'value']
    # Split 'Row' and 'Column' columns
    long_matrix[['chrom_x', 'start_x', 'end_x']] = long_matrix['Row'].str.split('[:-]', expand=True)
    long_matrix[['chrom_y', 'start_y', 'end_y']] = long_matrix['Column'].str.split('[:-]', expand=True)


    # Convert start1, end1, start2, end2 to numeric
    long_matrix[['start_x', 'end_x', 'start_y', 'end_y']] = long_matrix[['start_x', 'end_x', 'start_y', 'end_y']].apply(pd.to_numeric)

    # Reorder the columns
    long_matrix = long_matrix[['chrom_x', 'start_x', 'end_x', 'chrom_y', 'start_y', 'end_y', 'value']]

    return long_matrix

"""
specify directory path to save the output files
and create it if it does not exist
"""

directory_path = str(Path(args.output_dir) / '04_cool_files' / f'cool_files_npmi_permutation_test_{args.chr}') + '/'
os.makedirs(directory_path, exist_ok=True)
os.makedirs(str(Path(args.output_dir) / '04_cool_files'), exist_ok=True)

print("\n" + "="*40)


def save_long_matrix_nan_5(chr, resolution, data='npmi_matrix_hap1', data_matrix=npmi_matrix_hap1, directory_path=directory_path):
    npmi_matrix_df = add_chr_column(data_matrix, chr = chr, resolution = resolution)
    long_matrix = create_long_matrix(npmi_matrix_df)

    # generate filename for cool file for hap2
    output_file_path_long_mat = [f'{directory_path}permutation_test_results_{data}_{chr}_long_matrix_nan.tsv.gz']
    output_file_path_long_mat = ''.join(output_file_path_long_mat)

    with gzip.open(output_file_path_long_mat, 'wt', encoding='utf-8') as f:
        long_matrix.to_csv(f, sep='\t', index=False)

    print(f"Output {data} saved to {output_file_path_long_mat}")


    long_matrix['value'] = long_matrix['value'].fillna(5)

    output_file_path_nan = [f'{directory_path}permutation_test_results_{data}_{chr}_long_matrix_nan_is_five.tsv.gz']
    output_file_path_nan = ''.join(output_file_path_nan)
    with gzip.open(output_file_path_nan, 'wt', encoding='utf-8') as f:
        long_matrix.to_csv(f, sep='\t', index=False)

    print(f"Output {data} with NaN = 5 saved to {output_file_path_nan}")

# p_value_log_filter_npmi
save_long_matrix_nan_5(chr=args.chr, resolution = args.resolution, data='p_value_log_filter_npmi', data_matrix=p_value_log_filter_npmi, directory_path=directory_path)
# p_value_log_filter_thresh_npmi
save_long_matrix_nan_5(chr=args.chr, resolution = args.resolution, data='p_value_log_filter_thresh_npmi', data_matrix=p_value_log_filter_thresh_npmi, directory_path=directory_path)
# p_value_log_filter_thresh_npmi_directionality
save_long_matrix_nan_5(chr=args.chr, resolution = args.resolution, data='p_value_log_filter_thresh_npmi_directionality', data_matrix=p_value_log_filter_thresh_npmi_directionality, directory_path=directory_path)
# p_value_npmi_2sided_updated
save_long_matrix_nan_5(chr=args.chr, resolution = args.resolution, data='p_value_npmi_2sided_updated', data_matrix=p_value_npmi_2sided_updated, directory_path=directory_path)
# npmi_matrix_hap1
save_long_matrix_nan_5(chr=args.chr, resolution = args.resolution, data='npmi_matrix_hap1', data_matrix=npmi_matrix_hap1, directory_path=directory_path)
# npmi_matrix_hap2
save_long_matrix_nan_5(chr=args.chr, resolution = args.resolution, data='npmi_matrix_hap2', data_matrix=npmi_matrix_hap2, directory_path=directory_path)
# npmi_matrix_hap2
save_long_matrix_nan_5(chr=args.chr, resolution = args.resolution, data='npmi_matrix_both', data_matrix=npmi_matrix_both, directory_path=directory_path)
print("="*40 + "\n")
# #############


"""
create bed file with all the genomic windows (e.g. chr1 0 50000) 
"""
print("\n" + "="*40)

output_bed_file_path = str(Path(args.output_dir) / '04_cool_files' / f'res_{args.resolution}_genomic_windows.bed')

# File lock ensures only one parallel chromosome task generates the BED file;
# all others wait and then read the completed file.
lock_path = output_bed_file_path + '.lock'
with open(lock_path, 'w') as lock_file:
    fcntl.flock(lock_file, fcntl.LOCK_EX)
    if os.path.exists(output_bed_file_path):
        # Another task already wrote it while we were waiting for the lock
        bed_df = pd.read_csv(output_bed_file_path, sep='\t', header=None)
        print(f"Data loaded from {output_bed_file_path}")
    else:
        # First task to acquire the lock: download and write
        url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"

        response = requests.get(url)
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
                end = min(start + int(args.resolution), size)
                bed_data.append([chrom, start, end])
                start += int(args.resolution)

        bed_df = pd.DataFrame(bed_data, columns=['Chromosome', '0', str(args.resolution)])
        bed_df.to_csv(output_bed_file_path, sep='\t', index=False, header=False)
        print(f"Output saved to {output_bed_file_path}")
# lock released here; all waiting tasks now proceed to read the completed file
print("="*40 + "\n")


"""
generate cool format file
Adapted to keep NaN values and replace them with value 5

ADJUST TO EITHER GENERATE COOL FILES FOR HAP1, HAP2, BOTH OR FOR PERMUTATION TEST RESULTS

"""
#########################
# 
def create_cool_file(data='npmi_matrix_hap1'):
    cool_file_path = f'{directory_path}permutation_test_results_{data}_{args.chr}_nan_is_5.cool'
    output_file_path_nan = f'{directory_path}permutation_test_results_{data}_{args.chr}_long_matrix_nan_is_five.tsv.gz'

    # Bash command to create cool file
    cool_command = f"zcat {output_file_path_nan} | grep -v start_x | cooler load -f bg2 --count-as-float --assembly hg38 --input-copy-status duplex {shlex.quote(output_bed_file_path)} - {shlex.quote(cool_file_path)}"
    # Call the Bash command using subprocess
    subprocess.run(cool_command, check=True, shell=True, executable='/bin/bash')
    print(f"cool file has been generated for {data}: {cool_file_path}")    

# p_value_log_filter_npmi
create_cool_file(data='p_value_log_filter_npmi')
# p_value_log_filter_thresh_npmi
create_cool_file(data='p_value_log_filter_thresh_npmi')
# p_value_log_filter_thresh_npmi_directionality
create_cool_file(data='p_value_log_filter_thresh_npmi_directionality')
# p_value_npmi_2sided_updated
create_cool_file(data='p_value_npmi_2sided_updated')
# npmi_matrix_hap1
create_cool_file(data='npmi_matrix_hap1')
# npmi_matrix_hap2
create_cool_file(data='npmi_matrix_hap2')
# npmi_matrix_hap2
create_cool_file(data='npmi_matrix_both')
