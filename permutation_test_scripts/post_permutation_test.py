"""
Code by Claudia Robens
"""

import numpy as np
import pandas as pd
import pickle # used to save and load data
import scipy as sp
import glob # used to autocomplete filenames
import os
import copy
from scipy.ndimage import gaussian_filter


# find pkl file in directory based on the pattern "_{chrom}_" and load it using pickle load
def find_and_load_pkl(directory_path, chrom):
    # Define the pattern to search for .pkl files
    pattern = f"*_{chrom}_*.pkl"
    
    # Find the .pkl file in the subdirectories
    pkl_files = glob.glob(os.path.join(directory_path, pattern), recursive=True)
    
    # Check if any .pkl files were found
    if pkl_files:
        # Load the first .pkl file found
        print(pkl_files)
        with open(pkl_files[0], 'rb') as file:
            perm_results = pickle.load(file)
        return perm_results
    else:
        print(f"No .pkl files found with pattern {pattern} in {directory_path}")
        return None


def gaussian_filter_nan(arr, sigma=1, truncate=4.0, eps=1e-6):
    mask = np.isfinite(arr)
    arr_filled = np.where(mask, arr, 0.0)
    arr_filtered = sp.ndimage.gaussian_filter(arr_filled, sigma=sigma, truncate=truncate)
    mask_filtered = sp.ndimage.gaussian_filter(mask.astype(float), sigma=sigma, truncate=truncate)
    
    result = np.full_like(arr, np.nan)
    safe = (mask_filtered > eps) & np.isfinite(mask_filtered) & np.isfinite(arr_filtered)
    np.divide(arr_filtered, mask_filtered, out=result, where=safe)   # no invalid divide happens
    result[~mask] = np.nan # mask values that where originally nan
    result_sym = 0.5 * (result + result.T) #keep symmetry in the result
    np.fill_diagonal(result_sym, np.diag(result))  # preserve original diagonal
    return result_sym


def process_permutation_results(perm_results, chrom, bins_rm_f, gaussian_kernel_size=1 ):
    # Read the bins to be removed
    bins_rm = pd.read_csv(bins_rm_f, sep='\t', index_col = 0)
    bins_rm_filtered = bins_rm[bins_rm['chrom'] == chrom]
    tbr = perm_results['subset_segtable_hap1'].index.get_level_values('start').isin(bins_rm_filtered['start'])  # to be removed

    # CALCULATE UPDATED P-VALUE WITH CORRECTION
    perm_results['p_value_npmi_2sided_updated'] = np.minimum(1,2*np.minimum(((perm_results['perm_greater']+perm_results['perm_equal']+1 ) / (perm_results['perm_greater']+perm_results['perm_smaller']+perm_results['perm_equal']+1)), ((perm_results['perm_smaller']+perm_results['perm_equal']+1 ) / (perm_results['perm_greater']+perm_results['perm_smaller']+perm_results['perm_equal']+1))))

    # mask perm_results['p_value_npmi_2sided'] with nan mask of perm_results['npmi_diff_nan_mask']
    filtered_array = copy.deepcopy(perm_results)
    filtered_array['p_value_npmi_2sided_updated'][filtered_array['npmi_diff_nan_mask']] = np.nan
    
    
    # mask p-value matrix based on WDF curation
    filtered_array['p_value_npmi_2sided_updated'][:, np.array(tbr)] = np.nan  # set columns to NaN
    filtered_array['p_value_npmi_2sided_updated'][np.array(tbr), :] = np.nan 
    
    #log transformation of p-values
    filtered_array['p_value_log_npmi'] = -np.log10(filtered_array['p_value_npmi_2sided_updated'])

    # filtering using Gaussian filter
    filtered_array['p_value_log_filter_npmi'] = gaussian_filter_nan(filtered_array['p_value_log_npmi'], sigma = gaussian_kernel_size)
    
    # add directionality of the contacts by using 'perm_greater' and 'perm_smaller' stored in perm_results
    filtered_array['greater_mask'] = filtered_array['perm_greater'] < filtered_array['perm_smaller'] # stronger contact on hap1
    filtered_array['smaller_mask'] = filtered_array['perm_smaller'] < filtered_array['perm_greater'] # stronger contact on hap2
    
    filtered_array['directionality'] = filtered_array['npmi_hap1'] - filtered_array['npmi_hap2'] # stronger contact on hap2
    directionality_mask = filtered_array['directionality'] < 0

    filtered_array['p_value_log_filter_npmi_directionality'] = copy.deepcopy(filtered_array['p_value_log_filter_npmi'])
    filtered_array["p_value_log_filter_npmi_directionality"][directionality_mask] = -np.abs(filtered_array["p_value_log_filter_npmi_directionality"][directionality_mask]) # positive values stronger on hap1, negative stronger on hap2

    return filtered_array

def process_permutation_results_thresh(perm_results, chrom, thresh_table, bins_rm_f,subset_segtable_hap1, gaussian_kernel_size=1):
    # Read the bins to be removed
    bins_rm = pd.read_csv(bins_rm_f, sep='\t', index_col = 0)
    bins_rm_filtered = bins_rm[bins_rm['chrom'] == chrom]
    tbr = subset_segtable_hap1.index.get_level_values('start').isin(bins_rm_filtered['start'])  # to be removed # cannot use the subset_segtable_hap1 here, because it is only a matrix without indices in the perm_results

    # CALCULATE UPDATED P-VALUE WITH CORRECTION
    perm_results['p_value_npmi_2sided_updated'] = np.minimum(1,2*np.minimum(((perm_results['perm_greater']+perm_results['perm_equal']+1 ) / (perm_results['perm_greater']+perm_results['perm_smaller']+perm_results['perm_equal']+1)), ((perm_results['perm_smaller']+perm_results['perm_equal']+1 ) / (perm_results['perm_greater']+perm_results['perm_smaller']+perm_results['perm_equal']+1))))
    # mask perm_results['p_value_npmi_2sided'] with nan mask of perm_results['npmi_diff_nan_mask']
    filtered_array = copy.deepcopy(perm_results)
    filtered_array['p_value_npmi_2sided_updated'][filtered_array['npmi_diff_nan_mask']] = np.nan
    
    
    # mask p-value matrix based on WDF curation
    filtered_array['p_value_npmi_2sided_updated'][:, np.array(tbr)] = np.nan  # set columns to NaN
    filtered_array['p_value_npmi_2sided_updated'][np.array(tbr), :] = np.nan 
    
    #log transformation of p-values
    filtered_array['p_value_log_npmi'] = -np.log10(filtered_array['p_value_npmi_2sided_updated'])
    
    # smoothing using Gaussian filter
    filtered_array['p_value_log_filter_npmi'] = gaussian_filter_nan(filtered_array['p_value_log_npmi'], sigma = gaussian_kernel_size)
    
    # filtering based on thresholded
    filtered_array['p_value_log_filter_thresh_npmi'] = copy.deepcopy(filtered_array['p_value_log_filter_npmi'])
    
    # include thresholding step based on knee point from threshold table
    threshold_df = pd.read_csv(thresh_table, sep = "\t")
    threshold_df = threshold_df.set_index('chromosome')
    threshold = threshold_df.loc[chrom, 'knee']
    message = f"Threshold for {chrom}: {threshold}"
    print(message)
    filtered_array['p_value_log_filter_thresh_npmi'][filtered_array['p_value_log_filter_thresh_npmi'] < threshold] = 0


    # add directionality of the contacts by using 'perm_greater' and 'perm_smaller' stored in perm_results
    filtered_array['greater_mask'] = filtered_array['perm_greater'] < filtered_array['perm_smaller'] # stronger contact on hap1
    filtered_array['smaller_mask'] = filtered_array['perm_smaller'] < filtered_array['perm_greater'] # stronger contact on hap2
    
    filtered_array['directionality'] = filtered_array['npmi_hap1'] - filtered_array['npmi_hap2'] # stronger contact on hap2
    directionality_mask = filtered_array['directionality'] < 0

    # copy p_value_log_filter_thresh_npmi and create a new level were the directionality is stored 
    filtered_array['p_value_log_filter_thresh_npmi_directionality'] = copy.deepcopy(filtered_array['p_value_log_filter_thresh_npmi'])
    filtered_array["p_value_log_filter_thresh_npmi_directionality"][directionality_mask] = -np.abs(filtered_array["p_value_log_filter_thresh_npmi_directionality"][directionality_mask]) # positive values stronger on hap1, negative stronger on hap2

    return filtered_array