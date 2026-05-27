# `xgb_age_cage_gpu_full.py`

Functional README for the cage-aware mouse age prediction script.

This document describes the essential code behavior, required inputs, major processing steps, and the most important outputs written by the script.

## Script purpose

`xgb_age_cage_gpu_full.py` trains and evaluates an XGBoost regression model for mouse age prediction using abundance features while controlling leakage at the `Cage` level.

The script was designed to support:

- virus-only models
- bacteria-only models
- merged virus+bacteria models
- diet-specific or pooled `ALL` analyses

## Core design

The essential design principle is **cage-aware evaluation**.

This means:

- the held-out test set is split by `Cage`, not by individual sample
- outer CV is also group-aware at the `Cage` level
- no cage is shared between train and test within the holdout evaluation
- no cage is shared between train and validation within each CV fold

This reduces leakage from repeated or related samples housed in the same cage.

## What the script does

At a high level, the workflow is:

1. Read metadata and abundance table.
2. Align metadata and abundance by sample ID.
3. Optionally subset to a specific diet group, or use pooled `ALL`.
4. Construct age-based stratification labels.
   - for pooled `ALL`, the default stratification is `Diet + age stratum`
   - for diet-specific analyses, the default is age stratum only
5. Generate a **balanced cage-level holdout split**.
   - the split is chosen by randomized search over candidate cage sets
   - candidate splits are scored for balance in:
     - test fraction
     - mean age
     - strata composition
6. Optionally run **two-stage hyperparameter tuning** using:
   - `RandomizedSearchCV`
   - `StratifiedGroupKFold`
7. Run **outer CV** on the training set using cage-aware grouped folds.
8. Fit the final model on the full holdout-training set.
9. Predict the held-out test set.
10. Save split tables, prediction outputs, metrics, and a final XGBoost JSON model.

## Required inputs

The script requires:

- a metadata table
- an abundance table

### Required metadata columns

By default the script expects:

- `sample`
- `age.approx.wks`
- `mouse.ID`
- `Cage`
- `Diet`
- `Age bins`

These can be changed with CLI arguments such as:

- `--meta-sample-col`
- `--target-col`
- `--mouse-id-col`
- `--group-col`
- `--diet-col`
- `--age-strata-col`

### Abundance table requirements

The abundance table must contain one row per sample and numeric feature columns.

Supported sample-ID formats are:

- a column named `sample`
- an unnamed first column such as `Unnamed: 0`
- a user-specified column via `--abundance-sample-col`

All feature columns are forced to numeric. Non-numeric feature columns will trigger an error.

## Leakage-safe split logic

The holdout split is not a simple random sample split.

Instead, the script uses:

- `balanced_group_holdout_split(...)`

This function:

- treats `Cage` as the grouping variable
- randomly samples candidate test-cage sets
- scores each candidate using:
  - deviation from target test fraction
  - age distribution mismatch
  - strata imbalance
  - missing common strata in test

The best-scoring candidate becomes the final holdout split.

This is important because sklearn does not provide a built-in `StratifiedGroupShuffleSplit`.

## Cross-validation logic

Outer CV on the training set is performed using:

- `StratifiedGroupKFold` when possible
- `GroupKFold` fallback if stratified grouped CV fails

This preserves cage-level leakage control during model selection and internal evaluation.

## Hyperparameter tuning

Unless `--skip-greedy` is used, the script performs two-stage local tuning:

- **Stage 1**
  - randomized search around the base parameter set
- **Stage 2**
  - tighter randomized search around the best Stage 1 parameters

Tuning uses:

- `RandomizedSearchCV`
- grouped inner CV
- scoring = negative mean absolute error

## XGBoost / GPU behavior

The script uses XGBoost with:

- `tree_method="hist"`
- `device="cuda"` for GPU mode

Physical GPU selection is controlled by:

- `CUDA_VISIBLE_DEVICES`

through:

- `--gpu-id`

Example:

```bash
python xgb_age_cage_gpu_full.py \
  --meta do_meta.csv \
  --abundance MRGV_genus_all_features.csv \
  --out virus_age_ALL_cage_gpu.pkl \
  --diet-group ALL \
  --gpu-id 0 \
  --device cuda
```

For environment recreation, use:

- `xgb_gpu0.yml`

## Most important outputs

The script writes:

### 1. Result pickle

- `*.pkl`

This contains:

- final model
- feature names
- sample IDs
- training/test dataframes
- CV predictions and metrics
- holdout predictions and metrics

This is convenient, but it is not the most portable artifact.

### 2. Split tables directory

If `--out` is:

- `example.pkl`

the default split directory is:

- `example_split_tables/`

This directory contains the most portable reproducibility outputs:

- `holdout_train_clr.csv.gz`
- `holdout_test_clr.csv.gz`
- `holdout_train_meta.csv`
- `holdout_test_meta.csv`
- `holdout_train_samples.txt`
- `holdout_test_samples.txt`
- `holdout_test_predictions.csv`
- `holdout_test_metrics.json`
- `outer_cv_metrics.csv`
- `outer_cv_predictions.csv`
- `outer_cv_summary.json`
- `cv_splits.json`
- `final_xgb_model.json`
- `used_xgb_params.json`

### 3. Optional larger outputs

Depending on flags, the script may also save:

- per-fold train/validation CLR tables
- CV models inside the pickle
- tuning CV result tables inside the pickle

## Recommended reproducibility artifact

For public upload, the recommended artifact is:

- the split-tables directory
- plus `xgb_gpu0.yml`
- plus this script

This is preferable to uploading only the pickle because:

- the split tables are plain-text and portable
- the final model is stored as `final_xgb_model.json`
- the pickle is more sensitive to version mismatches in `pandas`, `xgboost`, and related packages

## Minimal command-line example

```bash
python xgb_age_cage_gpu_full.py \
  --meta do_meta.csv \
  --abundance MRGV_genus_abundance.csv \
  --out r358_virus_mouse_age_ALL_cage_gpu.pkl \
  --diet-group ALL \
  --gpu-id 0 \
  --device cuda \
  --outer-folds 10 \
  --inner-folds 5
```

## Notes for users

- The script assumes one row per sample in both metadata and abundance tables.
- Abundance features must already be numeric and preprocessed as intended by the analysis.
- The split package generated by this script is usually more suitable for GitHub, Zenodo, or supplementary upload than the pickle itself.

