"""
Code by Claudia Robens
"""

"""
load modules
"""
from generate_npmi import coordinates_from_location_string, get_region_from_location_string, calculate_NPMI, calculate_NPMI_steps, get_ticks, get_npmi_scale, read_segregation_table

import os
import pandas as pd
import numpy as np
import argparse # to parse command-line arguments
from datetime import datetime
import glob # used to autocomplete filenames
import pickle # used to save and load data
import logging
from multiprocessing import Pool # for parallel processing
from pathlib import Path


# define args
"""
Command-line arguments
"""
parser = argparse.ArgumentParser()
parser.add_argument('--output_dir', type=str, required=True, help='Base output directory; curated tables are read from output_dir/01_curated_segregation_tables/')
parser.add_argument('--cutoff', type=str, required=True, help='The cutoff value which was used for Co-Phasing pipeline, e.g. "10Mb" for 10 Mb cutoff.')
parser.add_argument('--resolution', type=int, required=True, help='The resolution which was used for Co-Phasing pipeline, e.g. "40000" for 40kb resolution.')
parser.add_argument('--num_perm', type=int, required=True, help='number of permutations ')
parser.add_argument('--chr', type=str, required=True, help='target chromosome chr1, chr2, etc.')
parser.add_argument('--pseudocount', type=int, default=0, help='pseudocount added to numerator and denominator when calculating NPMI, default is 0')
parser.add_argument("--out", required=True)  
parser.add_argument('--num_workers', type=int, default=10, help='number of parallel worker processes')


args = parser.parse_args()

print("\n" + "="*40)
start = datetime.now()
print(f"Start time: {start.strftime('%Y-%m-%d %H:%M:%S')}")

# Print arguments to stdout
print("\n" + "="*40)
print("Arguments applied:")
for arg, value in vars(args).items():
    print(f"{arg}: {value}")



"""
Read in segregation tables for Hap1, Hap2 haplotypes.
"""
print("="*40 + "\n")
print("Read in segregation tables for Hap1, Hap2 haplotypes")
# Hap1
curated_dir = str(Path(args.output_dir) / '01_curated_segregation_tables')
segregation_table_hap1 = read_segregation_table(Hap="hap1", file_path=curated_dir, resolution=args.resolution, cutoff=args.cutoff)
# Hap2
segregation_table_hap2 = read_segregation_table(Hap="hap2", file_path=curated_dir, resolution=args.resolution, cutoff=args.cutoff)

subset_segtable_hap1_df = get_region_from_location_string(segregation_table_hap1, args.chr)
subset_segtable_hap1 = subset_segtable_hap1_df.values
subset_segtable_hap2 = get_region_from_location_string(segregation_table_hap2, args.chr)
subset_segtable_hap2 = subset_segtable_hap2.values


"""
Multiprocessing
"""
# original data
npmi_hap1  = calculate_NPMI_steps(subset_segtable_hap1, pseudocount = args.pseudocount)
npmi_hap2  = calculate_NPMI_steps(subset_segtable_hap2, pseudocount = args.pseudocount)

print("\n" + "="*40)
print("Calculated NPMI for original data.")

npmi_diff = npmi_hap1 - npmi_hap2 # calculate difference between haplotypes for original data
npmi_diff_nan_mask = np.isnan(npmi_diff) # create nan mask to keep track of nan values in original difference matrix
# store results in sparse matrix to reduce memory
size = npmi_diff.shape[0]

#permutation test
subset_segtable_joint = np.concatenate((subset_segtable_hap1, subset_segtable_hap2), axis=1) # combine both haplotypes for permutation test

indices = list(range(subset_segtable_hap1.shape[1]*2)) # number of samples in the two haplotypes combined
# set seed for reproducibility
np.random.seed(42)
# create permutations of the indices which are used later on for the permutation test
permutations = [np.random.permutation(indices) for _ in range(args.num_perm)]
# blocksize = args.num_perm//10
# blocks = [permutations[i*blocksize:(i+1)*blocksize] for i in range(10)] # split permutations into 10 blocks of blocksize number of permutations each
# 
# blocks = [permutations[i::10] for i in range(min(10, args.num_perm))]
# use the number of workers to split the permutations into blocks for parallel processing
blocks = [permutations[i::args.num_workers] for i in range(min(args.num_workers, args.num_perm))]
total = sum(len(b) for b in blocks)
sizes = [len(b) for b in blocks]
print("\n" + "="*40)
print(f"{len(blocks)} blocks defined for permutation.")
print(f"num_perm={args.num_perm:4d}: {len(blocks)} blocks, sizes={sizes}, total={total}, ok={total==args.num_perm}")


"""
function that performs permutation test for each block of permutations
"""
def permutation_test_thread(perms):
    perm_greater = np.zeros((size, size), dtype=np.uint16)
    perm_smaller = np.zeros((size, size), dtype=np.uint16)
    perm_equal = np.zeros((size, size), dtype=np.uint16)
    for perm in perms:
        # get two permuted haplotypes
        permuted_hap1 = subset_segtable_joint[:,perm[:len(perm)//2]]
        permuted_hap2 = subset_segtable_joint[:,perm[len(perm)//2:]]
        # calculate NPMI for permuted haplotypes
        npmi_perm_hap1 = calculate_NPMI_steps(permuted_hap1, pseudocount = args.pseudocount)
        npmi_perm_hap2 = calculate_NPMI_steps(permuted_hap2, pseudocount = args.pseudocount)
        # calculate difference between permuted haplotypes
        npmi_perm_diff = npmi_perm_hap1 - npmi_perm_hap2
        # mask values that are nan in one the two difference matrices in both matrices (17.09.2025)
        # Create mask where either matrix has NaN
        nan_mask = np.isnan(npmi_diff) | np.isnan(npmi_perm_diff)
        # Apply mask to matrix generated from permuted haplotypes
        npmi_perm_diff[nan_mask] = np.nan
        # perm greater than original for two-sided test
        perm_greater_mask = npmi_perm_diff > npmi_diff
        perm_greater += perm_greater_mask
        # perm smaller than original for two-sided test
        perm_smaller_mask = npmi_perm_diff < npmi_diff
        perm_smaller += perm_smaller_mask
        # perm equal to original for two-sided test (17.09.2025)
        perm_equal_mask = npmi_perm_diff == npmi_diff
        perm_equal += perm_equal_mask
    return perm_greater, perm_smaller, perm_equal
    

print("\n" + "="*40)
print("permutation test function defined.")

if __name__ == '__main__':
    # Create a pool with 10 processes
    # with Pool(processes=10) as pool:
    # use the number of workers specified in the command-line arguments to create a pool of processes
    with Pool(processes=args.num_workers) as pool:
        # Map compute_pxy_parallel function across all matrices
        results = pool.map(permutation_test_thread, blocks)
        pool.close()
        pool.join()
        perm_greater, perm_smaller, perm_equal = zip(*results)
        
print("\n" + "="*40)
print("pool function finished.")

"""
sum up results from all blocks
"""
perm_greater = np.sum(perm_greater, axis=0)
perm_smaller = np.sum(perm_smaller, axis=0)
perm_equal = np.sum(perm_equal, axis=0)

print("\n" + "="*40)
print("permutation test outputs combined.")

# Store the results in the dictionary
perm_results = {
    'npmi_hap1': npmi_hap1,
    'npmi_hap2': npmi_hap2,
    'npmi_diff_nan_mask': npmi_diff_nan_mask,
    'subset_segtable_hap1' : subset_segtable_hap1_df,
    'perm_greater': perm_greater, 
    'perm_smaller' : perm_smaller,
    'perm_equal' : perm_equal
}


    
# 2) Ensure Nextflow-visible deterministic output IN THE WORK DIR:
out_path = Path(args.out)

# ensure directory exists
out_path.parent.mkdir(parents=True, exist_ok=True)

# atomic write
tmp = out_path.with_suffix(out_path.suffix + ".tmp")
with open(tmp, 'wb') as f:
    pickle.dump(perm_results, f)

os.replace(tmp, out_path)

print(f"perm_results saved in {out_path}")
logging.info(f"perm_results saved in {out_path}")


print("="*40 + "\n")
end = datetime.now()
duration = end - start
print(f"End time: {end.strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Total runtime: {duration}")