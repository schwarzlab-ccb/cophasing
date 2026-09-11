import pandas as pd
import numpy as np
import copy
import argparse
from pathlib import Path
from generate_npmi import read_segregation_table


def get_chromosomes_list(segregation_table):
    """
    Note: Code by Alexander Kukalev
    """
    all_chromosomes = list(set(segregation_table.index.get_level_values('chrom')))
    all_chromosomes.sort()
    karyotype_set = []
    for chromosome in all_chromosomes:
        if (chromosome.find('_') == -1) and chromosome != "chrM":
            karyotype_set.append(chromosome)
    return karyotype_set


def get_WDF(segregation_table):
    """
    Calculate window detection frequency (WDF)
    
    Note: 
    Code by Alexander Kukalev
    Winick-Ng W, Kukalev A, Harabula I, Zea-Redondo L, Szabó D, Meijer M, et al. Cell-type specialization is encoded by specific chromatin topologies. Nature. 2021;599: 684–691. doi:10.1038/s41586-021-04081-2
    https://github.com/pombo-lab/WinickNg_Kukalev_Harabula_Nature_2021/blob/main/code/ZScore_differential_pipeline_GAMbrain.py
    """
    if len(segregation_table.index.names) == 1:
        segregation_table.set_index(['chrom', 'start', 'stop'], inplace=True)
    number_of_samples = len(segregation_table.columns)
    windows_sum = segregation_table.sum(axis=1).to_numpy()
    frequencies_list = windows_sum / number_of_samples
    WDF_df = pd.DataFrame(index=segregation_table.index)
    WDF_df['WDF'] = frequencies_list
    return WDF_df


def curate_segtable_with_mean_WDF_exclude_zeros(segregation_table, high_factor=1.5, low_factor=1.5):
     """    
        Curate a segregation table by removing windows with outlier WDF values.
        
        Notes:
        Code adapted from Alexander Kukalev
                
        Adapted by Claudia Robens: 
        split factor into separate high_factor and low_factor to allow asymmetric thresholds; 
        added defensive copy of WDF_values before zero-masking; 
        updated DataFrame.append()to _append().
    """
    segtable = segregation_table.copy()
    chrom_list = get_chromosomes_list(segtable)
    WDF_df = get_WDF(segtable)
    WDF_df.reset_index(inplace=True)
    segtable.reset_index(inplace=True, drop=True)
    curated_segtable = segtable.copy()
    removed_regions_df = pd.DataFrame(columns=['chrom', 'start', 'stop', 'WDF'], dtype=object)
    for chromosome in chrom_list:
        subset_WDF = WDF_df.query('chrom == @chromosome')
        WDF_values = subset_WDF['WDF'].to_numpy()
        WDF_values_eval = WDF_values.copy()
        WDF_values_eval[WDF_values_eval == 0] = np.nan
        mean_WDF = np.nanmean(WDF_values_eval)
        sigma_WDF = np.nanstd(WDF_values_eval)
        high_threshold = mean_WDF + high_factor * sigma_WDF
        low_threshold = mean_WDF - low_factor * sigma_WDF
        curated_WDF_index = subset_WDF.query('WDF <= @low_threshold or WDF >= @high_threshold')
        curated_segtable.loc[curated_WDF_index.index] = 0
        removed_regions_df = removed_regions_df._append(curated_WDF_index, ignore_index=True)
    curated_segtable.index = segregation_table.index
    total_windows = len(curated_segtable)
    removed_windows = len(removed_regions_df)
    percent = (removed_windows / total_windows) * 100
    print('Percent of removed windows is -> ' + str(percent))
    return curated_segtable, removed_regions_df, percent


def zero_rows(df, coords_df):
    idx_main = df.reset_index().set_index(['chrom', 'start']).index
    idx_to_zero = coords_df.set_index(['chrom', 'start']).index
    mask = idx_main.isin(idx_to_zero)
    df.loc[mask, :] = 0
    return df


parser = argparse.ArgumentParser(description='Curate segregation tables by mean WDF filtering')
parser.add_argument('--input_dir', type=str, required=True,
                    help='Directory containing raw segregation tables from the CoPhasing pipeline')
parser.add_argument('--output_dir', type=str, required=True,
                    help='Base output directory; curated tables are written to output_dir/01_curated_segregation_tables/')
parser.add_argument('--cutoff', type=str, required=True,
                    help='Cutoff value used in the CoPhasing pipeline, e.g. "10Mb"')
parser.add_argument('--resolution', type=int, required=True,
                    help='Resolution used in the CoPhasing pipeline, e.g. 40000')
parser.add_argument('--high_factor', type=float, default=1.8,
                    help='Upper WDF threshold factor (mean + high_factor * sigma), default 1.8')
parser.add_argument('--low_factor', type=float, default=1.8,
                    help='Lower WDF threshold factor (mean - low_factor * sigma), default 1.8')
args = parser.parse_args()

Path(args.output_dir).mkdir(parents=True, exist_ok=True)
output_subdir = Path(args.output_dir) / '01_curated_segregation_tables'
output_subdir.mkdir(parents=True, exist_ok=True)

print("Reading segregation tables...")
seg_table_df_hap1 = read_segregation_table(Hap='hap1', file_path=args.input_dir,
                                            resolution=args.resolution, cutoff=args.cutoff)
seg_table_df_hap2 = read_segregation_table(Hap='hap2', file_path=args.input_dir,
                                            resolution=args.resolution, cutoff=args.cutoff)
seg_table_df_both = read_segregation_table(Hap='both', file_path=args.input_dir,
                                            resolution=args.resolution, cutoff=args.cutoff)

print("Curating hap1...")
_, removed_hap1, pct_hap1 = curate_segtable_with_mean_WDF_exclude_zeros(
    seg_table_df_hap1, args.high_factor, args.low_factor)
print("Curating hap2...")
_, removed_hap2, pct_hap2 = curate_segtable_with_mean_WDF_exclude_zeros(
    seg_table_df_hap2, args.high_factor, args.low_factor)

# Cross-apply curation: zero out bins flagged in either haplotype
seg_table_df_hap1_cur = copy.deepcopy(seg_table_df_hap1)
seg_table_df_hap1_cur = zero_rows(seg_table_df_hap1_cur, removed_hap1)
seg_table_df_hap1_cur = zero_rows(seg_table_df_hap1_cur, removed_hap2)

seg_table_df_hap2_cur = copy.deepcopy(seg_table_df_hap2)
seg_table_df_hap2_cur = zero_rows(seg_table_df_hap2_cur, removed_hap2)
seg_table_df_hap2_cur = zero_rows(seg_table_df_hap2_cur, removed_hap1)

seg_table_df_both_cur = copy.deepcopy(seg_table_df_both)
seg_table_df_both_cur = zero_rows(seg_table_df_both_cur, removed_hap1)
seg_table_df_both_cur = zero_rows(seg_table_df_both_cur, removed_hap2)

# Save curated tables — naming matches read_segregation_table glob pattern:
# *_{cutoff}_cutoff.{resolution}.{hap}.segregation*
hap1_out = output_subdir / f"curated_{args.cutoff}_cutoff.{args.resolution}.hap1.segregation.tsv"
hap2_out = output_subdir / f"curated_{args.cutoff}_cutoff.{args.resolution}.hap2.segregation.tsv"
both_out = output_subdir / f"curated_{args.cutoff}_cutoff.{args.resolution}.both.segregation.tsv"

seg_table_df_hap1_cur.to_csv(hap1_out, sep='\t', index=True, header=True)
seg_table_df_hap2_cur.to_csv(hap2_out, sep='\t', index=True, header=True)
seg_table_df_both_cur.to_csv(both_out, sep='\t', index=True, header=True)

# Save combined list of bins removed in either haplotype
bins_rm_hap1_hap2 = pd.merge(
    removed_hap1[['chrom', 'start', 'stop']],
    removed_hap2[['chrom', 'start', 'stop']],
    on=['chrom', 'start', 'stop'],
    how='outer'
)
bins_rm_hap1_hap2.to_csv(output_subdir / 'bins_rm_hap1_hap2.tsv', sep='\t', index=True, header=True)

print(f"Curated segregation tables saved to {output_subdir}")
print(f"  hap1: {hap1_out.name}")
print(f"  hap2: {hap2_out.name}")
print(f"  both: {both_out.name}")
print(f"  bins removed: bins_rm_hap1_hap2.tsv")
