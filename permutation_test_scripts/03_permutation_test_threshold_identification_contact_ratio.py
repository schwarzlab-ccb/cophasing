"""
Code by Claudia Robens
"""

from post_permutation_test import find_and_load_pkl, process_permutation_results

import os
import pandas as pd
import numpy as np
import argparse # to parse command-line arguments
import glob # used to autocomplete filenames
from kneed import KneeLocator # to find the knee point in a curve 
from pathlib import Path
import logging


"""
Command-line arguments
"""
parser = argparse.ArgumentParser()
parser.add_argument('--output_dir', type=str, required=True,
                    help='Base output directory containing all pipeline subdirectories')
parser.add_argument('--chromosomes', type=str, required=True,
                    help='Comma-separated list of chromosomes to process, e.g. "chr1,chr2,...,chr22"')
parser.add_argument('--gaussian_kernel_size', type=int, required=True,
                    help='Kernel size used for Gaussian filter')

args = parser.parse_args()
logging.info("args done.")

curated_dir = str(Path(args.output_dir) / '01_curated_segregation_tables')
permutation_test_dir = str(Path(args.output_dir) / '02_permutation_test')
bins_rm_file = str(Path(args.output_dir) / '01_curated_segregation_tables' / 'bins_rm_hap1_hap2.tsv')
threshold_output_dir = str(Path(args.output_dir) / '03_thresholds')

# Print arguments to stdout
print("\n" + "="*40)
print("Arguments applied:")
for arg, value in vars(args).items():
    print(f"{arg}: {value}")
print("="*40 + "\n")


"""
data
read in table with bins to be removed based on WDF curation
"""

chromosomes = args.chromosomes.split(',')

# Dictionary to store ratios for each chromosome
ratios_dict = {}
knee_points = {}

"""
CALCULATE RATIOS AND KNEE-POINT FOR ALL CHROMOSOMES
WITHOUT NAN VALUES AS TOTAL
"""
# Iterate over all chromosomes
for chrom in chromosomes:
    print(f"Processing {chrom}...")
    perm_results = find_and_load_pkl(directory_path=permutation_test_dir, chrom=chrom)

    filtered_array = process_permutation_results(perm_results=perm_results, chrom=chrom,
                                                  bins_rm_f=bins_rm_file, gaussian_kernel_size=args.gaussian_kernel_size)
    
    # Mask and process the array
    values = filtered_array['p_value_log_filter_npmi_directionality'].copy()
    lower_triangle_mask = np.tril(np.ones(values.shape), -1).astype(bool)
    filtered_array['p_value_log_filter_npmi_directionality_masked'] = np.where(
        lower_triangle_mask, np.nan, filtered_array['p_value_log_filter_npmi_directionality']
    )

    # Extract upper triangle values and remove NaNs
    upper_triangle_indices = np.triu_indices(filtered_array['p_value_log_filter_npmi_directionality_masked'].shape[0])
    p_value_log_filter_npmi_directionality_masked_values = filtered_array['p_value_log_filter_npmi_directionality_masked'][upper_triangle_indices]
    p_value_log_filter_npmi_directionality_masked_values_wo_nan = p_value_log_filter_npmi_directionality_masked_values[~np.isnan(p_value_log_filter_npmi_directionality_masked_values)]

    # Calculate thresholds and ratios
    threshold_step_size = 0.05
    thresholds = np.arange(0, np.nanmax(np.abs(filtered_array['p_value_log_filter_npmi_directionality'].flatten())), threshold_step_size)
    ratios = []
    for threshold in thresholds:
        # WITHOUT NAN AS TOTAL
        count_above_threshold = np.sum(np.abs(p_value_log_filter_npmi_directionality_masked_values_wo_nan) > threshold)
        total_non_nan = len(p_value_log_filter_npmi_directionality_masked_values_wo_nan)
        ratio_above_threshold = count_above_threshold / total_non_nan
        ratios.append(ratio_above_threshold)


    # Store the ratios for the current chromosome
    ratios_dict[chrom] = ratios
    if threshold_step_size == 0.1:
        # Use KneeLocator to find the elbow point
        # without nan values
        knee = KneeLocator(thresholds[1:], ratios[1:], curve='convex', direction='decreasing')
        knee_points[chrom] = knee.knee
    elif threshold_step_size == 0.05:
        knee = KneeLocator(thresholds[3:], ratios[3:], curve='convex', direction='decreasing')
        knee_points[chrom] = knee.knee
   

"""
WITHOUT NAN AS TOTAL
total number excluding masked nan values
"""
knee_df = pd.DataFrame(knee_points.items(), columns=['chromosome', 'knee'])
# Find the maximum length of the ratio lists
max_length = max(len(ratios) for ratios in ratios_dict.values())

# Extend each chromosome's ratio list with zeros to match the maximum length
for chrom in ratios_dict:
    ratios_dict[chrom] += [np.nan] * (max_length - len(ratios_dict[chrom]))

# Convert the dictionary to a DataFrame
ratios_df = pd.DataFrame(ratios_dict)

# add indices 0 to max_length*0.1
ratios_df.index = np.arange(0, max_length*threshold_step_size, threshold_step_size)

# Save the DataFrame to TSV files
threshold_ratio_total_wo_nan_output_dir = threshold_output_dir + '/'
if not os.path.exists(threshold_ratio_total_wo_nan_output_dir):
    os.makedirs(threshold_ratio_total_wo_nan_output_dir)

output_tsv_path = os.path.join(f"{threshold_ratio_total_wo_nan_output_dir}/ratios_table_threshold_steps.tsv")
ratios_df.to_csv(output_tsv_path, sep = '\t')


output_tsv_path = os.path.join(f"{threshold_ratio_total_wo_nan_output_dir}/knee_point_table_threshold_steps.tsv")
knee_df.to_csv(output_tsv_path, sep = '\t')

"""
WITHOUT NAN AS TOTAL
total number excluding masked nan values
"""
knee_df = knee_df.set_index('chromosome')

print("\n" + "="*40)
print("Thresholds calculated for all chromosomes using the knee point method.")
print("="*40 + "\n")


"""
calculate ratios for each chromosome separate for stronger on hap1 and hap2
"""
# get dataframe with knee points for each chromosome
threshold_df_file = os.path.join(threshold_ratio_total_wo_nan_output_dir, 'knee_point_table_threshold_steps.tsv')
threshold_df = pd.read_csv(threshold_df_file, sep='\t', index_col=0)
threshold_df = threshold_df.set_index('chromosome')

ratios_hap1_dict = {}
ratios_hap2_dict = {}
# Iterate over all chromosomes
for chrom in chromosomes:
    print(f"Processing {chrom}...")
    perm_results = find_and_load_pkl(directory_path=permutation_test_dir, chrom=chrom)
    filtered_array = process_permutation_results(perm_results=perm_results, chrom=chrom,
                                                  bins_rm_f=bins_rm_file, gaussian_kernel_size=args.gaussian_kernel_size)
    
    
    # Mask and process the array
    values = filtered_array['p_value_log_filter_npmi_directionality'].copy()
    lower_triangle_mask = np.tril(np.ones(values.shape), -1).astype(bool)
    filtered_array['p_value_log_filter_npmi_directionality_masked'] = np.where(
        lower_triangle_mask, np.nan, filtered_array['p_value_log_filter_npmi_directionality']
    )

    # Extract upper triangle values and remove NaNs
    upper_triangle_indices = np.triu_indices(filtered_array['p_value_log_filter_npmi_directionality_masked'].shape[0])
    p_value_log_filter_npmi_directionality_masked_values = filtered_array['p_value_log_filter_npmi_directionality_masked'][upper_triangle_indices]
    p_value_log_filter_npmi_directionality_masked_values_wo_nan = p_value_log_filter_npmi_directionality_masked_values[~np.isnan(p_value_log_filter_npmi_directionality_masked_values)]

    # get threshold value for respective chromosome
    threshold = threshold_df.loc[chrom, 'knee']
    # Calculate ratios
    # ALL WO NAN VALUES as TOTAL
    count_above_threshold = np.sum(p_value_log_filter_npmi_directionality_masked_values_wo_nan > threshold)
    total_non_nan = len(p_value_log_filter_npmi_directionality_masked_values_wo_nan)
    ratio_above_threshold = count_above_threshold / total_non_nan
    count_negative_threshold = np.sum(p_value_log_filter_npmi_directionality_masked_values_wo_nan < -threshold)
    ratio_negative_threshold = count_negative_threshold / total_non_nan
    # store results in the dictionary for each chromosome
    ratios_hap1_dict[chrom] = ratio_above_threshold
    ratios_hap2_dict[chrom] = ratio_negative_threshold
    
print("\n" + "="*40)
print("Ratios calculated for all chromosomes using the determined threshold per chromosome.")
print("="*40 + "\n")


"""save as tsv file"""
df_ratio = pd.DataFrame([
    {"chromosome": chr_, "hap1": ratios_hap1_dict.get(chr_, None), "hap2": ratios_hap2_dict.get(chr_, None)}
    for chr_ in chromosomes
])
df_threshold = pd.read_csv(f'{threshold_ratio_total_wo_nan_output_dir}/knee_point_table_threshold_steps.tsv', sep='\t', index_col=0)
df_threshold = df_threshold.rename(columns={"knee": "threshold"})
df_ratio_threshold = df_ratio.merge(df_threshold, on="chromosome", how="left")
df_ratio_threshold['combined_ratio'] = df_ratio_threshold['hap1'] + df_ratio_threshold['hap2']
df_ratio_threshold['hap1_percentage'] = (df_ratio_threshold['hap1'] * 100).round(2)
df_ratio_threshold['hap2_percentage'] = (df_ratio_threshold['hap2'] * 100).round(2)
df_ratio_threshold['combined_percentage'] = (df_ratio_threshold['hap1_percentage'] + df_ratio_threshold['hap2_percentage']).round(2)

df_ratio_threshold.to_csv(f'{threshold_ratio_total_wo_nan_output_dir}/ratio_differential_contacts_info.tsv', sep="\t", index=False)


