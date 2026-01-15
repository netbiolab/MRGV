#!/usr/bin/env python3

import os, sys, pickle, subprocess
import openpyxl
from multiprocessing import Process
import pandas as pd
from skbio.stats.composition import clr
from skbio import DistanceMatrix
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa
from matplotlib.patches import Ellipse
from collections import Counter
from skbio.stats.distance import permanova
from scipy import stats
from sklearn.metrics import mean_absolute_error
from xgboost import XGBRegressor
from typing import Tuple, List, Dict, Any
from sklearn.model_selection import (
    train_test_split, 
    RandomizedSearchCV,
    StratifiedGroupKFold,  
    GroupShuffleSplit      
)
from skbio.diversity import alpha_diversity
from scipy.spatial.distance import pdist, squareform
from itertools import combinations_with_replacement
import shap
import warnings
from scipy.stats import kendalltau
warnings.filterwarnings("ignore")


def _make_age_bins(y, n_bins: int = 10) -> np.ndarray:
    """
    Create quantile-based bins for continuous age 
    """
    y = np.asarray(y)
    percentiles = np.linspace(0, 100, n_bins + 1)
    bin_edges = np.unique(np.percentile(y, percentiles))
    if len(bin_edges) <= 2:
        return np.zeros_like(y, dtype=int)
    y_binned = np.digitize(y, bin_edges[1:-1], right=True)
    return y_binned.astype(int)


def _local_grid(center, deltas, vmin, vmax, as_int=False):
    """
    Build small local grid of candidate values around a center, clipped to [vmin, vmax].
    """
    vals = []
    for d in deltas:
        v = center + d
        if v < vmin or v > vmax:
            continue
        if as_int:
            v = int(round(v))
        vals.append(v)
    vals = sorted(set(vals))
    return vals


def _build_stage1_space(base_params: Dict[str, Any], n_features: int) -> Dict[str, List[Any]]:
    """
    First-round parameter search space:
    local neighborhood around base_params, aware of high dimensionality.
    """
    md = base_params.get("max_depth", 3)
    mcw = base_params.get("min_child_weight", 3)
    gamma = base_params.get("gamma", 1.0)
    subsample = base_params.get("subsample", 0.6)
    colsample = base_params.get("colsample_bytree", 0.3)
    reg_alpha = base_params.get("reg_alpha", 2.0)
    reg_lambda = base_params.get("reg_lambda", 2.0)
    lr = base_params.get("learning_rate", 0.01)
    n_estimators = base_params.get("n_estimators", 1200)

    # For very high dimensional data, keep colsample fairly small
    if n_features > 2000:
        col_max = 0.4
        col_min = 0.1
    else:
        col_max = 0.8
        col_min = 0.1

    param_distributions = {
        "max_depth": _local_grid(md, [-1, 0, 1], 2, 6, as_int=True),
        "min_child_weight": _local_grid(mcw, [-1, 0, 1, 2], 1, 10, as_int=True),
        "gamma": _local_grid(gamma, [-0.5, 0.0, 0.5], 0.0, 2.0, as_int=False),
        "subsample": _local_grid(subsample, [-0.1, 0.0, 0.1], 0.4, 1.0, as_int=False),
        "colsample_bytree": _local_grid(colsample, [-0.1, 0.0, 0.1], col_min, col_max, as_int=False),
        "reg_alpha": _local_grid(reg_alpha, [-1.0, 0.0, 1.0], 0.0, 5.0, as_int=False),
        "reg_lambda": _local_grid(reg_lambda, [-1.0, 0.0, 2.0], 0.0, 10.0, as_int=False),
        "learning_rate": _local_grid(lr, [-0.005, 0.0, 0.005], 0.005, 0.05, as_int=False),
        "n_estimators": _local_grid(n_estimators, [-400, 0, 400], 200, 2000, as_int=True),
    }
    return param_distributions


def _build_stage2_space(best_stage1: Dict[str, Any], n_features: int) -> Dict[str, List[Any]]:
    """
    Second-round (fine-tuning) search space:
    even tighter neighborhood around stage1 best.
    """
    md = best_stage1["max_depth"]
    mcw = best_stage1["min_child_weight"]
    gamma = best_stage1["gamma"]
    subsample = best_stage1["subsample"]
    colsample = best_stage1["colsample_bytree"]
    reg_alpha = best_stage1["reg_alpha"]
    reg_lambda = best_stage1["reg_lambda"]
    lr = best_stage1["learning_rate"]
    n_estimators = best_stage1["n_estimators"]

    if n_features > 2000:
        col_max = 0.4
        col_min = 0.1
    else:
        col_max = 0.8
        col_min = 0.1

    param_distributions = {
        "max_depth": _local_grid(md, [-1, 0, 1], 2, 6, as_int=True),
        "min_child_weight": _local_grid(mcw, [-1, 0, 1], 1, 10, as_int=True),
        "gamma": _local_grid(gamma, [-0.3, 0.0, 0.3], 0.0, 2.0, as_int=False),
        "subsample": _local_grid(subsample, [-0.05, 0.0, 0.05], 0.4, 1.0, as_int=False),
        "colsample_bytree": _local_grid(colsample, [-0.05, 0.0, 0.05], col_min, col_max, as_int=False),
        "reg_alpha": _local_grid(reg_alpha, [-0.5, 0.0, 0.5], 0.0, 5.0, as_int=False),
        "reg_lambda": _local_grid(reg_lambda, [-0.5, 0.0, 0.5], 0.0, 10.0, as_int=False),
        "learning_rate": _local_grid(lr, [-0.003, 0.0, 0.003], 0.003, 0.05, as_int=False),
        "n_estimators": _local_grid(n_estimators, [-200, 0, 200], 200, 2500, as_int=True),
    }
    return param_distributions
    
def run_xgboost_age_two_stage_tuning_5cv_10cv(
    abdf: pd.DataFrame,
    meta: pd.DataFrame,
    com_group: str,
    diet_group: str,
    mouse_id_col: str,  # <--- NEW ARGUMENT: Name of Mouse ID column
    base_best_params: Dict[str, Any],
    test_size: float = 0.2,
    random_state: int = 315,
    n_iter_stage1: int = 20,
    n_iter_stage2: int = 10,
    tune_fraction: float = 0.7,
    max_cpus: int = 100,
    skip_greedy: bool = False,   
) -> Dict[str, Any]:
    """
    Longitudinal-aware Two-stage XGBoost tuning.
    Uses GroupShuffleSplit for Train/Test and StratifiedGroupKFold for CV
    to prevent data leakage between timepoints of the same mouse.
    """

    # ---- build design matrix (filter by diet group) ----
    # Ensure MouseID is included in the concat so we can extract it later
    cols_to_keep = [com_group, "Diet", mouse_id_col]
    # Check if mouse_id_col exists in meta
    if mouse_id_col not in meta.columns:
        raise ValueError(f"Mouse ID column '{mouse_id_col}' not found in metadata!")

    mdf = pd.concat([abdf, meta[cols_to_keep]], axis=1)
    
    if diet_group != "ALL":
        mdf = mdf.loc[mdf["Diet"] == diet_group]
    assert len(mdf) != 0, f"No samples found for diet group '{diet_group}'"
    
    # Extract GROUPS (Mouse IDs)
    groups = mdf[mouse_id_col].values
    
    # Remove metadata cols to get just features (X) and target (y)
    mdf_features = mdf.drop(cols_to_keep, axis=1)
    feature_names = mdf_features.columns.tolist()
    
    X = mdf_features.values
    y = mdf[com_group].values

    n_samples, n_features = X.shape
    n_mice = len(np.unique(groups))
    print(f"Table shape for diet '{diet_group}': {n_samples} samples x {n_features} features")
    print(f"Total Unique Mice: {n_mice}")

    # ---- 1) Group-aware Train/Test split ----
    print("Performing Group-wise Train/Test split (longitudinal safe)...")
    # GroupShuffleSplit splits based on groups, not rows
    gss = GroupShuffleSplit(n_splits=1, test_size=test_size, random_state=random_state)
    train_idx, test_idx = next(gss.split(X, y, groups))
    
    X_train, X_test = X[train_idx], X[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]
    groups_train, groups_test = groups[train_idx], groups[test_idx]
    
    print(f"Train: {len(X_train)} samples ({len(np.unique(groups_train))} mice)")
    print(f"Test:  {len(X_test)} samples ({len(np.unique(groups_test))} mice)")

    # Create age bins for stratification usage in CV
    train_age_bins = _make_age_bins(y_train, n_bins=10)

    # ---- 1.5) Subsample training data for tuning (Group-aware) ----
    if 0 < tune_fraction < 1.0:
        gss_tune = GroupShuffleSplit(n_splits=1, train_size=tune_fraction, random_state=random_state)
        tune_idx, _ = next(gss_tune.split(X_train, y_train, groups_train))
        
        X_train_tune = X_train[tune_idx]
        y_train_tune = y_train[tune_idx]
        groups_tune = groups_train[tune_idx]
        bins_tune = train_age_bins[tune_idx]
    else:
        X_train_tune, y_train_tune = X_train, y_train
        groups_tune = groups_train
        bins_tune = train_age_bins

    print(f"Tuning on {X_train_tune.shape[0]} samples ({len(np.unique(groups_tune))} mice)")

    # base model (single-threaded to avoid nested parallelism explosion)
    base_model = XGBRegressor(
        objective="reg:squarederror",
        n_jobs=1,
        random_state=random_state,
        tree_method="hist",
    )

    # Placeholders for stage results
    best_params_stage1 = None
    best_params_stage2 = None
    best_cv_mae_stage1 = None
    best_cv_mae_stage2 = None

    # =====================================================
    #   OPTION A: Skip greedy search, use base_best_params
    # =====================================================
    if skip_greedy:
        print("Skipping greedy search. Using base_best_params directly for 10-fold CV.")
        model_params_for_cv = dict(base_best_params)

    # =====================================================
    #   OPTION B: Run Stage 1 + Stage 2 greedy search
    # =====================================================
    else:
        # ---- PREPARE GROUP-AWARE CV SPLITTER ----
        # StratifiedGroupKFold tries to balance 'bins_tune' while keeping 'groups_tune' intact
        skf_inner = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=random_state)
        
        # NOTE: We must generate the splits manually or pass groups to fit(). 
        # RandomizedSearchCV.fit(X, y, groups=groups) handles this correctly if cv=splitter.
        
        # ---- Stage 1: local search around global best ----
        param_stage1 = _build_stage1_space(base_best_params, n_features)
        print("Stage 1: RandomizedSearchCV (Group-Aware 5-fold)...")
        
        rnd_search1 = RandomizedSearchCV(
            estimator=base_model,
            param_distributions=param_stage1,
            n_iter=n_iter_stage1,
            scoring="neg_mean_absolute_error",
            n_jobs=max_cpus,
            cv=skf_inner,  # Pass the splitter object
            verbose=1,
            random_state=random_state,
            refit=True,
        )
        # IMPORTANT: Pass groups here!
        rnd_search1.fit(X_train_tune, y_train_tune, groups=groups_tune)
        
        best_params_stage1 = rnd_search1.best_params_
        best_cv_mae_stage1 = -rnd_search1.best_score_
        print("Best Stage 1 params:")
        for k, v in best_params_stage1.items():
            print(f"  {k}: {v}")
        print(f"Best Stage 1 inner-CV MAE: {best_cv_mae_stage1:.4f}")

        # ---- Stage 2: finer search around Stage 1 best ----
        param_stage2 = _build_stage2_space(best_params_stage1, n_features)
        print("Stage 2: fine-tuning around Stage 1 best params (5-fold)...")
        
        rnd_search2 = RandomizedSearchCV(
            estimator=base_model,
            param_distributions=param_stage2,
            n_iter=n_iter_stage2,
            scoring="neg_mean_absolute_error",
            n_jobs=max_cpus,
            cv=skf_inner, # Reuse splitter
            verbose=1,
            random_state=random_state,
            refit=True,
        )
        # IMPORTANT: Pass groups here!
        rnd_search2.fit(X_train_tune, y_train_tune, groups=groups_tune)
        
        best_params_stage2 = rnd_search2.best_params_
        best_cv_mae_stage2 = -rnd_search2.best_score_

        print("Best Stage 2 params:")
        for k, v in best_params_stage2.items():
            print(f"  {k}: {v}")
        print(f"Best Stage 2 inner-CV MAE: {best_cv_mae_stage2:.4f}")

        model_params_for_cv = dict(best_params_stage2)

    # ---- 2) Outer Group-Aware CV on full train ----
    print("Outer 10-fold Stratified Group CV with selected hyperparameters...")
    skf_outer = StratifiedGroupKFold(n_splits=10, shuffle=True, random_state=random_state)

    cv_mae_scores: List[float] = []
    
    # Manually loop to print fold results
    # split(X, y, groups)
    for n, (train_idx, val_idx) in enumerate(skf_outer.split(X_train, train_age_bins, groups_train)):
        print(f" CV fold {n + 1}/10...")
        X_tr, X_val = X_train[train_idx], X_train[val_idx]
        y_tr, y_val = y_train[train_idx], y_train[val_idx]

        model = XGBRegressor(
            objective="reg:squarederror",
            n_jobs=max_cpus,  # single-threaded model
            random_state=random_state,
            tree_method="hist",
            **model_params_for_cv,
        )
        model.fit(X_tr, y_tr)
        y_val_pred = model.predict(X_val)
        fold_mae = mean_absolute_error(y_val, y_val_pred)
        cv_mae_scores.append(fold_mae)
        print("  Fold %s MAE: %.4f" % (n + 1, fold_mae))

    cv_mae_mean = float(np.mean(cv_mae_scores))
    cv_mae_std = float(np.std(cv_mae_scores))

    # ---- 3) Final model on full training set ----
    final_model = XGBRegressor(
        objective="reg:squarederror",
        n_jobs=max_cpus,
        random_state=random_state,
        tree_method="hist",
        **model_params_for_cv,
    )
    final_model.fit(X_train, y_train)

    # ---- 4) Evaluate on held-out test set (New Mice) ----
    y_pred_test = final_model.predict(X_test)
    test_mae = mean_absolute_error(y_test, y_pred_test)
    print("Outer CV Mean MAE: %.4f ± %.4f" % (cv_mae_mean, cv_mae_std))
    print("Test MAE (Held-out Mice): %.4f" % test_mae)

    return {
        "final_model": final_model,
        "cv_scores": cv_mae_scores,
        "cv_mean": cv_mae_mean,
        "cv_std": cv_mae_std,
        "test_mae": test_mae,
        "y_test": y_test,
        "y_pred_test": y_pred_test,
        "X_train": X_train,
        "X_test": X_test,
        "feature_names": feature_names,
        "best_params_stage1": best_params_stage1,
        "best_params_stage2": best_params_stage2,
        "stage1_cv_mae": best_cv_mae_stage1,
        "stage2_cv_mae": best_cv_mae_stage2,
        "ab_df": abdf,
        "meta_df": meta,
        "com_group": com_group,
        "diet_group": diet_group,
        "used_params_for_cv": model_params_for_cv,
        "skip_greedy": skip_greedy,
    }

if __name__ == "__main__":
    base_best_params = {
        "subsample":        0.6,
        "reg_lambda":       2.0,
        "reg_alpha":        2.0,
        "n_estimators":     1200,
        "min_child_weight": 3,
        "max_depth":        3,
        "learning_rate":    0.01,
        "gamma":            1.0,
        "colsample_bytree": 0.3,
    }

    meta = pd.read_csv(sys.argv[1], sep="\t", index_col = 0)
    abdf = pd.read_csv(sys.argv[2], sep="\t", index_col = 0)
    svpath = sys.argv[3]
    diet_group = sys.argv[4]  # e.g., "40"
    dtype = sys.argv[5]
    threads = int(sys.argv[6])
    greedy = bool(sys.argv[7])
    
    # --- Check for Mouse ID argument, default to "MouseID" if not provided ---
    try:
        mouse_id_col = sys.argv[8]
    except IndexError:
        # CHANGE THIS if your column name is different
        mouse_id_col = "MouseID"
        print(f"WARNING: No Mouse ID column provided. Defaulting to '{mouse_id_col}'.")

    print("------------------------------------")
    print(f"Running diet group: {diet_group}, data type: {dtype}")
    print(f"Using Mouse ID column: {mouse_id_col}")
    print("------------------------------------")
    
    res = run_xgboost_age_two_stage_tuning_5cv_10cv(
        abdf=abdf,
        meta=meta,
        com_group="age.approx.wks",
        diet_group=diet_group,
        mouse_id_col="mouse.ID",  # <--- Pass it here
        base_best_params=base_best_params,
        n_iter_stage1=50,
        n_iter_stage2=30,
        tune_fraction=1.0,
        max_cpus=threads,
        test_size=0.2,
        skip_greedy=greedy,
    )
    pickle.dump(res, open(svpath, "wb"))
