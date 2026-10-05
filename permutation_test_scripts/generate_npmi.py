import pandas as pd
import numpy as np
import glob

def coordinates_from_location_string (region):
    """
    Note: Code by Alexander Kukalev
    """
    chrom = region.split(':')[0]
    if len (region.split(':')) == 1:
        start, stop = '', ''
    else:
        coordinates = region.split(':')[1]
        start = coordinates.split('-') [0]
        stop = coordinates.split ('-') [1]
        start = int(start.replace(',', ''))
        stop = int(stop.replace(',', ''))
    if start > stop:
        raise Exception ('Check coordinates of your region! It can not end before the start')
    return chrom, start, stop

def get_region_from_location_string (segregation_table, region):
    """
    Note: Code by Alexander Kukalev

    Adapted during review: bins are selected by actual overlap with the requested region,
    and the range check compares the requested stop against the chromosome end instead of
    the start of its last bin — the latter rejected every region reaching the final,
    truncated bin of a chromosome.
    """
    if len (segregation_table.index.names) == 1: # Check if supplied segregation table has multiindex
        segregation_table.set_index(['chrom','start','stop'], inplace = True)
    chrom, start, stop = coordinates_from_location_string (region)
    seg_chrome = segregation_table.loc [chrom] # Subset for chromosome
    if start=='' and stop =='':
        return seg_chrome
    start_values = np.array (seg_chrome.index.get_level_values('start'))
    stop_values = np.array (seg_chrome.index.get_level_values('stop'))
    if stop > np.max (stop_values):
        raise Exception ('Check coordinates of your region! The chromosome is shorter than specified range')
    overlapping = (stop_values > start) & (start_values < stop) # Every bin touching the region
    if not overlapping.any():
        raise Exception ('Check coordinates of your region! No bins overlap the requested range')
    return seg_chrome [overlapping]

def calculate_NPMI (segregation_table, region):
    """    
        Notes:
        Adapted from:
        Winick-Ng W, Kukalev A, Harabula I, Zea-Redondo L, Szabó D, Meijer M, et al. Cell-type specialization is encoded by specific chromatin topologies. Nature. 2021;599: 684–691. doi:10.1038/s41586-021-04081-2
        https://github.com/pombo-lab/WinickNg_Kukalev_Harabula_Nature_2021/blob/main/code/ZScore_differential_pipeline_GAMbrain.py
                
        Adapted by Claudia Robens: 
        replaced len(seg_matrix[0]) with seg_matrix.shape[1], and dot() with @ operator.
    """
    subset_segtable = get_region_from_location_string (segregation_table, region)
    seg_matrix = subset_segtable.values
    M = seg_matrix.shape[1]
    pxy = (seg_matrix @ seg_matrix.T) / M #calculate p(x,y) matrix i.e coseg matrix divided my M
    #p(x)p(y) matrix  - i.e detection frequency with the dot product of it's transposed self equals an N*N matrix
    pxpy = seg_matrix.sum(1).reshape(-1,1)/M * seg_matrix.sum(1)/M
    PMI = np.log2(pxy/pxpy) #define PMI matrix
    NPMI = PMI/-np.log2(pxy) #bound values between -1 and 1
    return NPMI

# with pseudocount added to numerator and denominator
def calculate_NPMI_steps (seg_matrix, pseudocount=1):
    """    
        Notes:
        Adapted from:
        Winick-Ng W, Kukalev A, Harabula I, Zea-Redondo L, Szabó D, Meijer M, et al. Cell-type specialization is encoded by specific chromatin topologies. Nature. 2021;599: 684–691. doi:10.1038/s41586-021-04081-2
        https://github.com/pombo-lab/WinickNg_Kukalev_Harabula_Nature_2021/blob/main/code/ZScore_differential_pipeline_GAMbrain.py
                
        Adapted by Claudia Robens: 
        added pseudocount parameter applied to numerator and denominator of pxy and pxpy,
        replaced len(seg_matrix[0]) with seg_matrix.shape[1], and dot() with @ operator.
    """
    M = seg_matrix.shape[1]
    # calculate p(x,y) matrix i.e coseg matrix divided my M
    pxy = ((seg_matrix @ seg_matrix.T)+pseudocount) / (M +pseudocount)
    #p(x)p(y) matrix  - i.e detection frequency with the dot product of it's transposed self equals an N*N matrix
    pxpy = (seg_matrix.sum(1)+pseudocount).reshape(-1,1)/(M+pseudocount) * (seg_matrix.sum(1)+pseudocount)/(M+pseudocount) 
    PMI = np.log2(pxy/pxpy) #define PMI matrix
    NPMI = PMI/-np.log2(pxy) #bound values between -1 and 1
    return NPMI


def get_ticks (segregation_table, region):
    """
    Code by Alexander Kukalev
    """
    chrom, start, stop = coordinates_from_location_string (region)
    seg_chrome = segregation_table.loc [chrom] # Subset for chromosome
    if start=='' and stop =='':
        start_values = np.array (seg_chrome.index.get_level_values('start'))
        start = 0
        stop = np.max (start_values)
        total_bins = len (seg_chrome)
    else:
        seg_table = get_region_from_location_string (segregation_table, region)
        total_bins = len (seg_table)
    middle = int(start + ((stop - start) / 2))
    middle_bin = int (total_bins / 2)
    return middle_bin, total_bins, start, middle, stop

def get_npmi_scale (matrix):
    """
    Note: Code by Alexander Kukalev
    """
    matrix_con = np.concatenate (matrix)
    matrix_remove_nans = matrix_con[~np.isnan(matrix_con)]
    minscale = 0
    maxscale = np.percentile (matrix_remove_nans, 99)
    return minscale, maxscale

def read_segregation_table(Hap, file_path, cutoff, resolution):
    """
    Reads the segregation table for the specified haplotype and returns it as a pandas DataFrame.
    Raises FileNotFoundError if no table matches, ValueError if several do.
    Args:
        Hap: 'hap1', 'hap2', or 'both'. Defines the haplotype for which the segregation table should be read.
        file_path: The path to the directory containing the segregation table files.
        cutoff: The CoPhasing distance threshold value used in the Co-Phasing pipeline and chosen for analysis(e.g., '10Mb').
        resolution: The resolution used in the Co-Phasing pipeline and chosen for analysis (e.g., 40000 for 40kb).
    Returns:
        segregation_table: pandas DataFrame with the segregation table.

    """
    # Prefer tables whose name carries the cutoff; fall back to any table at this
    # resolution so output of a run that does not encode the cutoff still works.
    pattern = f"{file_path}/*_{str(cutoff)}_cutoff.{str(resolution)}.{str(Hap)}.segregation*"
    matching_files = sorted(glob.glob(pattern))
    if not matching_files:
        pattern = f"{file_path}/*.{str(resolution)}.{str(Hap)}.segregation*"
        matching_files = sorted(glob.glob(pattern))
    if len(matching_files) > 1:
        raise ValueError(
            f"{Hap}: {len(matching_files)} segregation tables match '{pattern}', "
            f"cannot choose between them: {matching_files}"
        )
    if not matching_files:
        raise FileNotFoundError(f"{Hap}: no segregation table matches '{pattern}'")
    hap = matching_files[0]
    print(f"{str(Hap)} segregation table: '{hap}'")
    segregation_table = pd.read_csv(hap, sep='\t', index_col=[0, 1, 2])
    return segregation_table