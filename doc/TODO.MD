# Open items — permutation test pipeline

Raised during the review of `feature/permutation-test` (2026-10-05) and deliberately
left unfixed: each one either changes published numbers or is a redesign. Evidence
below comes from runs against `test_data/` at resolution 1000.

## 1. The p-value denominator varies per cell

`post_permutation_test._p_values_and_smoothing` divides by `perm_greater +
perm_smaller + perm_equal + 1`, not by `num_perm + 1`. Permutations that produce
NaN are counted in none of the three, because NaN comparisons are always false
(`02_permutation_test_chromosome_level.py`, where `npmi_perm_diff[nan_mask] = np.nan`).

Measured on chr5 with `num_perm=50`: the valid count ranged from 24 to 48 across
cells, and the resulting p-value differed from the standard estimator by up to 0.42
in 33 of the 40 unmasked cells.

Two consequences:

* The attainable p-value floor is `2/(n+1)`, which ranged from 0.041 to 0.059 on that
  chromosome. Step 03 picks one knee threshold per chromosome and applies it to
  `-log10(p)` values that are therefore not on a common scale.
* The dropped permutations are not missing at random — they are the degenerate splits
  that put all detections of a bin into one haplotype — so the retained null is a
  biased subsample of the intended one.

`n` is never recorded, so none of this is visible downstream.

**Options:** record the per-cell valid count in the pkl and carry it into the threshold
table (no published number changes); or divide by `num_perm + 1`, treating NaN
permutations as non-extreme (every existing p-value changes).

## 2. `low_factor` also decides whether empty bins are dropped

`01_segregation_table_curation_by_mean_WDF.py` excludes zeros when computing the mean
and sigma, but still thresholds against them. All-zero bins are removed only while
`mean - low_factor * sigma >= 0`.

Measured on chr5: at `low_factor=1.8` the 22 empty bins are removed; at `3.0` the
threshold goes negative and all 22 survive. A chromosome that is entirely empty gives
`nanmean -> NaN`, so nothing is removed there at all.

One parameter therefore governs two unrelated decisions, and the behaviour flips per
chromosome depending on the spread of the data.

**Options:** drop all-zero bins unconditionally and apply the WDF outlier thresholds
separately (changes which bins are curated wherever the threshold currently goes
negative); or keep the coupling and document it.

## 3. Pipeline output is written out-of-band

All four processes write their real results into `params.output_dir` from inside the
Python scripts. Only the step 2 pkl is a declared Nextflow output; everything else
reaches `output_dir` as a side effect, and the declared outputs are log files.

Verified with a stub workflow of the same topology: on `-resume` all tasks cache on
work-directory state alone. Wiping or editing `output_dir` and resuming therefore
reports success while producing nothing, or leaves stale results in place.

Related: step 03 reads the pkls from the published `02_permutation_test` directory
rather than the files Nextflow stages into its work directory, so the dependency is
expressed only as a cache-busting token and not as real dataflow. Two concurrent runs
sharing an `output_dir` will also overwrite each other.

**Options:** have each step declare its outputs and consume staged inputs, so `-resume`
and concurrent runs become correct; or document that `output_dir` must not be touched
between runs and that `-resume` is unsafe if it is.

## 4. hg38 chrom.sizes is downloaded at runtime

`04_generate_cool_file_permutation_test_results.py` fetches
`https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes` while the
pipeline runs. There is no offline fallback and no way to point at a local copy, so the
pipeline needs outbound HTTP from every compute node. The `fcntl.flock` guard that
serialises the download across parallel chromosome tasks is unreliable over NFS.

The assembly is hardcoded in three places: the URL, the `autosomes = chr1..chr22` list,
and `cooler load --assembly hg38`.

**Options:** add a `--chrom_sizes` parameter and fall back to the download when it is
not given; or require a local file and drop the download entirely.

## Also noted, not raised as decisions

* The pipeline declares no `memory` or `time` directives and no executor profile. For
  chr1 at 40 kb (~6200 bins) each worker holds three uint32 count matrices of ~154 MB
  plus float64 NPMI matrices of ~308 MB; at `num_workers = 10` that reaches tens of GB.
* `pcp-env.yaml` pins no versions, so a fresh `conda env create` can resolve a different
  toolchain than the one the code was written against.
* Step 02 relies on `fork` to inherit module-level globals into the worker pool, so it
  is Linux-only.
