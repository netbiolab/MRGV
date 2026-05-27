#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cage-aware mouse age prediction with XGBoost GPU.

Main purpose
------------
1. Align metadata and CLR abundance table by sample ID.
2. Perform leakage-safe holdout split using Cage as the group variable.
3. Balance holdout split approximately by age strata, and by Diet+age strata for pooled ALL analysis.
4. Run optional two-stage hyperparameter tuning with StratifiedGroupKFold.
5. Run outer CV on the training set with Cage-aware StratifiedGroupKFold.
6. Fit final model on full training set and evaluate held-out Cage-level test set.
7. Save holdout train/test CLR tables, metadata, sample IDs, CV split metadata, and predictions.

GPU behavior
------------
Use CUDA_VISIBLE_DEVICES to choose the physical GPU, while passing only device='cuda'
to XGBoost. Example:

    python xgb_age_cage_gpu_full.py \
      --meta do_meta.csv \
      --abundance "MRGV_genus_all_features_ab_df(1).csv" \
      --out mouse_age_ALL_cage_gpu.pkl \
      --diet-group ALL \
      --gpu-id 3 \
      --device cuda

Internally:
    os.environ["CUDA_VISIBLE_DEVICES"] = "3"
    XGBRegressor(device="cuda", tree_method="hist")

Author: ChatGPT-assisted script for cage-aware longitudinal mouse age prediction.
"""

from __future__ import annotations

import argparse
import json
import os
import pickle
import sys
import warnings
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold, RandomizedSearchCV, StratifiedGroupKFold
from xgboost import XGBRegressor

warnings.filterwarnings("ignore")


# -----------------------------------------------------------------------------
# I/O helpers
# -----------------------------------------------------------------------------

def read_table_auto(path: str | Path) -> pd.DataFrame:
    """Read CSV/TSV using pandas delimiter inference."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    # sep=None with python engine handles comma/tab fairly robustly.
    return pd.read_csv(path, sep=None, engine="python")


def prepare_metadata(
    meta_path: str | Path,
    sample_col: str = "sample",
) -> pd.DataFrame:
    """Read metadata and set sample column as index."""
    meta = read_table_auto(meta_path)

    if sample_col not in meta.columns:
        raise ValueError(
            f"sample_col='{sample_col}' not found in metadata. "
            f"Available columns: {list(meta.columns)}"
        )

    if meta[sample_col].duplicated().any():
        dup = meta.loc[meta[sample_col].duplicated(), sample_col].head(10).tolist()
        raise ValueError(f"Duplicated sample IDs in metadata. Examples: {dup}")

    meta = meta.copy()
    meta[sample_col] = meta[sample_col].astype(str)
    meta = meta.set_index(sample_col, drop=False)
    return meta


def prepare_abundance(
    abundance_path: str | Path,
    sample_col: Optional[str] = None,
) -> pd.DataFrame:
    """
    Read abundance table and set sample IDs as index.

    Supported formats:
      1. A column named 'sample'.
      2. A first unnamed column, e.g. 'Unnamed: 0', containing sample IDs.
      3. User-provided sample_col.
    """
    ab = read_table_auto(abundance_path)

    if sample_col is not None:
        if sample_col not in ab.columns:
            raise ValueError(
                f"abundance sample_col='{sample_col}' not found. "
                f"Available first columns: {list(ab.columns[:10])}"
            )
        index_col = sample_col
    elif "sample" in ab.columns:
        index_col = "sample"
    elif len(ab.columns) > 0 and str(ab.columns[0]).startswith("Unnamed"):
        index_col = ab.columns[0]
    else:
        # Fallback: assume first column is sample ID if non-numeric/object-like.
        first_col = ab.columns[0]
        if ab[first_col].dtype == object:
            index_col = first_col
        else:
            raise ValueError(
                "Could not infer sample ID column in abundance table. "
                "Provide --abundance-sample-col explicitly."
            )

    if ab[index_col].duplicated().any():
        dup = ab.loc[ab[index_col].duplicated(), index_col].head(10).tolist()
        raise ValueError(f"Duplicated sample IDs in abundance table. Examples: {dup}")

    ab = ab.copy()
    ab[index_col] = ab[index_col].astype(str)
    ab = ab.set_index(index_col, drop=True)

    # Force features to numeric. Fail early if there are unexpected non-numeric columns.
    non_numeric = []
    for c in ab.columns:
        if not pd.api.types.is_numeric_dtype(ab[c]):
            try:
                ab[c] = pd.to_numeric(ab[c])
            except Exception:
                non_numeric.append(c)
    if non_numeric:
        raise ValueError(
            f"Non-numeric abundance feature columns detected: {non_numeric[:20]}"
        )

    return ab


def align_meta_abundance(
    meta: pd.DataFrame,
    ab: pd.DataFrame,
    target_col: str,
    group_col: str,
    mouse_id_col: str,
    diet_col: str,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Inner-align metadata and abundance by sample ID index, then validate columns."""
    required_cols = [target_col, group_col, mouse_id_col, diet_col]
    missing = [c for c in required_cols if c not in meta.columns]
    if missing:
        raise ValueError(
            f"Missing required metadata columns: {missing}. "
            f"Available columns: {list(meta.columns)}"
        )

    common = meta.index.intersection(ab.index)
    missing_in_ab = meta.index.difference(ab.index)
    missing_in_meta = ab.index.difference(meta.index)

    print("[Alignment]")
    print(f"  Metadata samples:  {meta.shape[0]:,}")
    print(f"  Abundance samples: {ab.shape[0]:,}")
    print(f"  Common samples:    {len(common):,}")
    print(f"  Missing in abundance: {len(missing_in_ab):,}")
    print(f"  Missing in metadata:  {len(missing_in_meta):,}")

    if len(common) == 0:
        raise ValueError("No overlapping sample IDs between metadata and abundance table.")

    # Keep metadata order after filtering to common samples.
    meta2 = meta.loc[common].copy()
    ab2 = ab.loc[common].copy()

    # Drop samples with missing target or group/diet fields.
    before = meta2.shape[0]
    keep_mask = meta2[required_cols].notna().all(axis=1)
    meta2 = meta2.loc[keep_mask].copy()
    ab2 = ab2.loc[meta2.index].copy()
    after = meta2.shape[0]
    if before != after:
        print(f"  Dropped samples with missing required metadata: {before - after:,}")

    meta2[target_col] = pd.to_numeric(meta2[target_col], errors="coerce")
    keep_mask = meta2[target_col].notna()
    if not keep_mask.all():
        print(f"  Dropped samples with non-numeric target: {(~keep_mask).sum():,}")
        meta2 = meta2.loc[keep_mask].copy()
        ab2 = ab2.loc[meta2.index].copy()

    return meta2, ab2


# -----------------------------------------------------------------------------
# Split helpers
# -----------------------------------------------------------------------------

def make_quantile_bins(y: Sequence[float], n_bins: int = 10) -> np.ndarray:
    """Create quantile-based bins for continuous target values."""
    y = np.asarray(y, dtype=float)
    if len(np.unique(y)) <= 1:
        return np.zeros_like(y, dtype=int)

    edges = np.unique(np.nanpercentile(y, np.linspace(0, 100, n_bins + 1)))
    if len(edges) <= 2:
        return np.zeros_like(y, dtype=int)

    return np.digitize(y, edges[1:-1], right=True).astype(int)


def make_strata(
    meta: pd.DataFrame,
    target_col: str,
    diet_col: str,
    age_strata_col: Optional[str],
    diet_group: str,
    n_age_bins: int,
    force_diet_in_strata: bool,
) -> pd.Series:
    """
    Build stratification labels.

    For ALL pooled data, default is Diet + age stratum.
    For diet-specific data, default is age stratum only.
    """
    if age_strata_col is not None and age_strata_col in meta.columns:
        age_strata = meta[age_strata_col].astype(str).fillna("NA")
    else:
        bins = make_quantile_bins(meta[target_col].values, n_bins=n_age_bins)
        age_strata = pd.Series(bins.astype(str), index=meta.index, name="quantile_age_bin")

    use_diet = force_diet_in_strata or (str(diet_group).upper() == "ALL")
    if use_diet:
        strata = meta[diet_col].astype(str) + "__" + age_strata.astype(str)
    else:
        strata = age_strata.astype(str)

    strata.name = "strata"
    return strata


def summarize_split_distribution(
    y: np.ndarray,
    strata: np.ndarray,
    groups: np.ndarray,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
) -> Dict[str, Any]:
    """Summarize train/test split distribution."""
    def vc(labels: np.ndarray) -> Dict[str, int]:
        s = pd.Series(labels.astype(str))
        return {str(k): int(v) for k, v in s.value_counts().sort_index().items()}

    train_groups = set(groups[train_idx].astype(str))
    test_groups = set(groups[test_idx].astype(str))

    out = {
        "n_train_samples": int(len(train_idx)),
        "n_test_samples": int(len(test_idx)),
        "n_train_groups": int(len(train_groups)),
        "n_test_groups": int(len(test_groups)),
        "n_shared_groups": int(len(train_groups.intersection(test_groups))),
        "train_mean_age": float(np.mean(y[train_idx])),
        "test_mean_age": float(np.mean(y[test_idx])),
        "train_median_age": float(np.median(y[train_idx])),
        "test_median_age": float(np.median(y[test_idx])),
        "train_min_age": float(np.min(y[train_idx])),
        "train_max_age": float(np.max(y[train_idx])),
        "test_min_age": float(np.min(y[test_idx])),
        "test_max_age": float(np.max(y[test_idx])),
        "train_strata_counts": vc(strata[train_idx]),
        "test_strata_counts": vc(strata[test_idx]),
    }
    return out


def _split_score(
    y: np.ndarray,
    strata: np.ndarray,
    test_idx: np.ndarray,
    target_test_size: float,
    all_strata_counts: pd.Series,
) -> float:
    """
    Score candidate holdout split.

    Lower is better. Penalizes:
      - sample test fraction deviation
      - mean age deviation
      - strata fraction deviation
      - missing common strata in test
    """
    n = len(y)
    test_mask = np.zeros(n, dtype=bool)
    test_mask[test_idx] = True

    obs_frac = float(test_mask.mean())
    sample_frac_penalty = abs(obs_frac - target_test_size) * 10.0

    y_all_mean = float(np.mean(y))
    y_all_sd = float(np.std(y)) if float(np.std(y)) > 0 else 1.0
    age_penalty = abs(float(np.mean(y[test_mask])) - y_all_mean) / y_all_sd

    test_counts = pd.Series(strata[test_mask].astype(str)).value_counts()
    strata_penalty = 0.0
    missing_penalty = 0.0

    for s, total_count in all_strata_counts.items():
        total_count = int(total_count)
        if total_count == 0:
            continue
        observed = int(test_counts.get(s, 0))
        desired = total_count * target_test_size
        # Normalize by sqrt(total_count) to avoid rare strata completely dominating.
        strata_penalty += ((observed - desired) ** 2) / max(total_count, 1)
        # Penalize absence in test for reasonably represented strata.
        if total_count >= max(5, int(round(1.0 / max(target_test_size, 1e-6)))) and observed == 0:
            missing_penalty += 2.0

    return float(sample_frac_penalty + age_penalty + strata_penalty + missing_penalty)


def balanced_group_holdout_split(
    y: np.ndarray,
    strata: np.ndarray,
    groups: np.ndarray,
    test_size: float = 0.2,
    random_state: int = 358,
    n_candidates: int = 5000,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    """
    Select a Cage/group-level holdout split by randomized search.

    This is used because sklearn does not provide StratifiedGroupShuffleSplit.
    The function samples many group subsets and chooses the one with best age/strata balance.
    """
    rng = np.random.default_rng(random_state)

    y = np.asarray(y, dtype=float)
    strata = np.asarray(strata).astype(str)
    groups = np.asarray(groups).astype(str)

    unique_groups = np.array(sorted(pd.unique(groups)))
    n_groups = len(unique_groups)
    if n_groups < 2:
        raise ValueError("Need at least 2 unique groups for holdout split.")

    n_test_groups = int(round(n_groups * test_size))
    n_test_groups = max(1, min(n_test_groups, n_groups - 1))

    group_to_indices = {g: np.where(groups == g)[0] for g in unique_groups}
    all_strata_counts = pd.Series(strata).value_counts()

    best_score = np.inf
    best_test_idx: Optional[np.ndarray] = None
    best_test_groups: Optional[np.ndarray] = None

    # Include deterministic offset candidates plus randomized candidates.
    candidates_checked = 0
    for i in range(max(1, n_candidates)):
        if i == 0:
            permuted = unique_groups.copy()
        else:
            permuted = rng.permutation(unique_groups)

        test_groups = np.sort(permuted[:n_test_groups])
        test_idx = np.concatenate([group_to_indices[g] for g in test_groups])
        score = _split_score(y, strata, test_idx, test_size, all_strata_counts)
        candidates_checked += 1

        if score < best_score:
            best_score = score
            best_test_idx = np.sort(test_idx)
            best_test_groups = test_groups.copy()

    assert best_test_idx is not None
    test_idx = best_test_idx
    train_idx = np.setdiff1d(np.arange(len(y)), test_idx)

    report = summarize_split_distribution(y, strata, groups, train_idx, test_idx)
    report.update(
        {
            "split_method": "balanced_group_holdout_random_search",
            "target_test_size": float(test_size),
            "n_candidate_splits_checked": int(candidates_checked),
            "best_score": float(best_score),
            "test_groups": [str(x) for x in best_test_groups],
        }
    )
    return train_idx, test_idx, report


def make_cv_splits(
    X: np.ndarray,
    strata: np.ndarray,
    groups: np.ndarray,
    n_splits: int,
    random_state: int,
) -> Tuple[List[Tuple[np.ndarray, np.ndarray]], str]:
    """
    Build StratifiedGroupKFold splits. Fall back to GroupKFold if needed.
    """
    n_unique_groups = len(pd.unique(groups))
    if n_unique_groups < 2:
        raise ValueError("Need at least 2 unique groups for CV.")

    n_splits = min(n_splits, n_unique_groups)
    if n_splits < 2:
        raise ValueError("n_splits became < 2 after considering unique groups.")

    try:
        cv = StratifiedGroupKFold(
            n_splits=n_splits,
            shuffle=True,
            random_state=random_state,
        )
        splits = list(cv.split(X, strata, groups))
        return splits, "StratifiedGroupKFold"
    except Exception as e:
        print(f"[WARNING] StratifiedGroupKFold failed: {e}")
        print("[WARNING] Falling back to GroupKFold. Age/Diet balance may be weaker.")
        cv = GroupKFold(n_splits=n_splits)
        splits = list(cv.split(X, groups=groups))
        return splits, "GroupKFold_fallback"


# -----------------------------------------------------------------------------
# Model helpers
# -----------------------------------------------------------------------------

def local_grid(center: Any, deltas: Sequence[Any], vmin: float, vmax: float, as_int: bool = False) -> List[Any]:
    vals = []
    for d in deltas:
        v = center + d
        if v < vmin or v > vmax:
            continue
        if as_int:
            v = int(round(v))
        vals.append(v)
    return sorted(set(vals))


def build_stage1_space(base_params: Dict[str, Any], n_features: int) -> Dict[str, List[Any]]:
    """First-round local search space around base parameters."""
    md = base_params.get("max_depth", 3)
    mcw = base_params.get("min_child_weight", 3)
    gamma = base_params.get("gamma", 1.0)
    subsample = base_params.get("subsample", 0.6)
    colsample = base_params.get("colsample_bytree", 0.3)
    reg_alpha = base_params.get("reg_alpha", 2.0)
    reg_lambda = base_params.get("reg_lambda", 2.0)
    lr = base_params.get("learning_rate", 0.01)
    n_estimators = base_params.get("n_estimators", 1200)

    if n_features > 2000:
        col_min, col_max = 0.1, 0.4
    else:
        col_min, col_max = 0.1, 0.8

    return {
        "max_depth": local_grid(md, [-1, 0, 1], 2, 6, as_int=True),
        "min_child_weight": local_grid(mcw, [-1, 0, 1, 2], 1, 10, as_int=True),
        "gamma": local_grid(gamma, [-0.5, 0.0, 0.5], 0.0, 2.0),
        "subsample": local_grid(subsample, [-0.1, 0.0, 0.1], 0.4, 1.0),
        "colsample_bytree": local_grid(colsample, [-0.1, 0.0, 0.1], col_min, col_max),
        "reg_alpha": local_grid(reg_alpha, [-1.0, 0.0, 1.0], 0.0, 5.0),
        "reg_lambda": local_grid(reg_lambda, [-1.0, 0.0, 2.0], 0.0, 10.0),
        "learning_rate": local_grid(lr, [-0.005, 0.0, 0.005], 0.003, 0.05),
        "n_estimators": local_grid(n_estimators, [-400, 0, 400], 200, 2500, as_int=True),
    }


def build_stage2_space(best_params: Dict[str, Any], n_features: int) -> Dict[str, List[Any]]:
    """Second-round tighter search space around stage1 best parameters."""
    md = best_params["max_depth"]
    mcw = best_params["min_child_weight"]
    gamma = best_params["gamma"]
    subsample = best_params["subsample"]
    colsample = best_params["colsample_bytree"]
    reg_alpha = best_params["reg_alpha"]
    reg_lambda = best_params["reg_lambda"]
    lr = best_params["learning_rate"]
    n_estimators = best_params["n_estimators"]

    if n_features > 2000:
        col_min, col_max = 0.1, 0.4
    else:
        col_min, col_max = 0.1, 0.8

    return {
        "max_depth": local_grid(md, [-1, 0, 1], 2, 6, as_int=True),
        "min_child_weight": local_grid(mcw, [-1, 0, 1], 1, 10, as_int=True),
        "gamma": local_grid(gamma, [-0.3, 0.0, 0.3], 0.0, 2.0),
        "subsample": local_grid(subsample, [-0.05, 0.0, 0.05], 0.4, 1.0),
        "colsample_bytree": local_grid(colsample, [-0.05, 0.0, 0.05], col_min, col_max),
        "reg_alpha": local_grid(reg_alpha, [-0.5, 0.0, 0.5], 0.0, 5.0),
        "reg_lambda": local_grid(reg_lambda, [-0.5, 0.0, 0.5], 0.0, 10.0),
        "learning_rate": local_grid(lr, [-0.003, 0.0, 0.003], 0.003, 0.05),
        "n_estimators": local_grid(n_estimators, [-200, 0, 200], 200, 3000, as_int=True),
    }


def make_xgb_regressor(
    params: Dict[str, Any],
    random_state: int,
    threads: int,
    device: str,
) -> XGBRegressor:
    """Create an XGBRegressor with requested CPU/GPU device."""
    model_kwargs = {
        "objective": "reg:squarederror",
        "random_state": random_state,
        "n_jobs": threads,
        "tree_method": "hist",
        **params,
    }

    if device.lower() != "cpu":
        # User requested: pass device='cuda', not cuda:<ordinal>.
        model_kwargs["device"] = device
    else:
        model_kwargs["device"] = "cpu"

    return XGBRegressor(**model_kwargs)


def regression_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    """Compute common regression metrics."""
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    mse = float(mean_squared_error(y_true, y_pred))
    rmse = float(np.sqrt(mse))
    mae = float(mean_absolute_error(y_true, y_pred))
    r2 = float(r2_score(y_true, y_pred)) if len(np.unique(y_true)) > 1 else float("nan")
    corr = float(np.corrcoef(y_true, y_pred)[0, 1]) if len(y_true) > 1 else float("nan")
    return {"mae": mae, "rmse": rmse, "r2": r2, "pearson_r": corr}


def two_stage_tuning(
    X_train_tune: np.ndarray,
    y_train_tune: np.ndarray,
    strata_tune: np.ndarray,
    groups_tune: np.ndarray,
    base_params: Dict[str, Any],
    n_features: int,
    random_state: int,
    threads: int,
    device: str,
    inner_folds: int,
    n_iter_stage1: int,
    n_iter_stage2: int,
    search_n_jobs: int,
) -> Dict[str, Any]:
    """Run two-stage RandomizedSearchCV with Cage-aware CV splits."""
    base_model = make_xgb_regressor(
        params={},
        random_state=random_state,
        threads=threads,
        device=device,
    )

    inner_splits, inner_cv_name = make_cv_splits(
        X=X_train_tune,
        strata=strata_tune,
        groups=groups_tune,
        n_splits=inner_folds,
        random_state=random_state,
    )

    print(f"[Tuning] Inner CV: {inner_cv_name}, folds={len(inner_splits)}")

    stage1_space = build_stage1_space(base_params, n_features=n_features)
    print("[Tuning] Stage 1 RandomizedSearchCV")
    search1 = RandomizedSearchCV(
        estimator=base_model,
        param_distributions=stage1_space,
        n_iter=n_iter_stage1,
        scoring="neg_mean_absolute_error",
        n_jobs=search_n_jobs,
        cv=inner_splits,
        verbose=1,
        random_state=random_state,
        refit=True,
        error_score="raise",
    )
    search1.fit(X_train_tune, y_train_tune)
    best_stage1 = dict(search1.best_params_)
    best_stage1_mae = float(-search1.best_score_)

    print("[Tuning] Best Stage 1 params:")
    for k, v in best_stage1.items():
        print(f"  {k}: {v}")
    print(f"[Tuning] Best Stage 1 MAE: {best_stage1_mae:.4f}")

    stage2_space = build_stage2_space(best_stage1, n_features=n_features)
    print("[Tuning] Stage 2 RandomizedSearchCV")
    search2 = RandomizedSearchCV(
        estimator=base_model,
        param_distributions=stage2_space,
        n_iter=n_iter_stage2,
        scoring="neg_mean_absolute_error",
        n_jobs=search_n_jobs,
        cv=inner_splits,
        verbose=1,
        random_state=random_state,
        refit=True,
        error_score="raise",
    )
    search2.fit(X_train_tune, y_train_tune)
    best_stage2 = dict(search2.best_params_)
    best_stage2_mae = float(-search2.best_score_)

    print("[Tuning] Best Stage 2 params:")
    for k, v in best_stage2.items():
        print(f"  {k}: {v}")
    print(f"[Tuning] Best Stage 2 MAE: {best_stage2_mae:.4f}")

    return {
        "best_params": best_stage2,
        "best_params_stage1": best_stage1,
        "best_params_stage2": best_stage2,
        "stage1_cv_mae": best_stage1_mae,
        "stage2_cv_mae": best_stage2_mae,
        "inner_cv_name": inner_cv_name,
        "inner_folds": len(inner_splits),
        "stage1_cv_results": pd.DataFrame(search1.cv_results_),
        "stage2_cv_results": pd.DataFrame(search2.cv_results_),
    }


# -----------------------------------------------------------------------------
# Save helpers
# -----------------------------------------------------------------------------

def ensure_dir(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_series_txt(values: Sequence[str], path: str | Path) -> None:
    with open(path, "w") as f:
        for v in values:
            f.write(f"{v}\n")


def save_json(obj: Dict[str, Any], path: str | Path) -> None:
    with open(path, "w") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)


def save_holdout_tables(
    outdir: Path,
    X_df: pd.DataFrame,
    meta: pd.DataFrame,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
) -> Dict[str, str]:
    """Save holdout train/test CLR tables and metadata."""
    train_samples = X_df.index[train_idx].astype(str).tolist()
    test_samples = X_df.index[test_idx].astype(str).tolist()

    paths = {
        "holdout_train_clr": str(outdir / "holdout_train_clr.csv.gz"),
        "holdout_test_clr": str(outdir / "holdout_test_clr.csv.gz"),
        "holdout_train_meta": str(outdir / "holdout_train_meta.csv"),
        "holdout_test_meta": str(outdir / "holdout_test_meta.csv"),
        "holdout_train_samples": str(outdir / "holdout_train_samples.txt"),
        "holdout_test_samples": str(outdir / "holdout_test_samples.txt"),
    }

    X_df.iloc[train_idx].to_csv(paths["holdout_train_clr"], compression="gzip")
    X_df.iloc[test_idx].to_csv(paths["holdout_test_clr"], compression="gzip")
    meta.iloc[train_idx].to_csv(paths["holdout_train_meta"])
    meta.iloc[test_idx].to_csv(paths["holdout_test_meta"])
    write_series_txt(train_samples, paths["holdout_train_samples"])
    write_series_txt(test_samples, paths["holdout_test_samples"])

    return paths


def save_cv_fold_info(
    outdir: Path,
    X_train_df: pd.DataFrame,
    meta_train: pd.DataFrame,
    splits: List[Tuple[np.ndarray, np.ndarray]],
    groups_train: np.ndarray,
    strata_train: np.ndarray,
    save_cv_clr_tables: bool,
) -> List[Dict[str, Any]]:
    """Save CV fold sample IDs and metadata. Optionally save CLR tables."""
    cv_records: List[Dict[str, Any]] = []

    for fold_id, (tr_idx, val_idx) in enumerate(splits, start=1):
        fold_dir = ensure_dir(outdir / f"cv_fold_{fold_id:02d}")

        train_samples = X_train_df.index[tr_idx].astype(str).tolist()
        val_samples = X_train_df.index[val_idx].astype(str).tolist()

        write_series_txt(train_samples, fold_dir / "train_samples.txt")
        write_series_txt(val_samples, fold_dir / "val_samples.txt")
        meta_train.iloc[tr_idx].to_csv(fold_dir / "train_meta.csv")
        meta_train.iloc[val_idx].to_csv(fold_dir / "val_meta.csv")

        clr_paths = {}
        if save_cv_clr_tables:
            train_clr_path = fold_dir / "train_clr.csv.gz"
            val_clr_path = fold_dir / "val_clr.csv.gz"
            X_train_df.iloc[tr_idx].to_csv(train_clr_path, compression="gzip")
            X_train_df.iloc[val_idx].to_csv(val_clr_path, compression="gzip")
            clr_paths = {
                "train_clr": str(train_clr_path),
                "val_clr": str(val_clr_path),
            }

        rec = {
            "fold": int(fold_id),
            "n_train_samples": int(len(tr_idx)),
            "n_val_samples": int(len(val_idx)),
            "n_train_groups": int(len(set(groups_train[tr_idx].astype(str)))),
            "n_val_groups": int(len(set(groups_train[val_idx].astype(str)))),
            "n_shared_groups": int(
                len(set(groups_train[tr_idx].astype(str)).intersection(set(groups_train[val_idx].astype(str))))
            ),
            "train_samples_path": str(fold_dir / "train_samples.txt"),
            "val_samples_path": str(fold_dir / "val_samples.txt"),
            "train_meta_path": str(fold_dir / "train_meta.csv"),
            "val_meta_path": str(fold_dir / "val_meta.csv"),
            "train_strata_counts": pd.Series(strata_train[tr_idx].astype(str)).value_counts().sort_index().to_dict(),
            "val_strata_counts": pd.Series(strata_train[val_idx].astype(str)).value_counts().sort_index().to_dict(),
            **clr_paths,
        }
        cv_records.append(rec)

    save_json({"folds": cv_records}, outdir / "cv_splits.json")
    return cv_records


# -----------------------------------------------------------------------------
# Main workflow
# -----------------------------------------------------------------------------

def run_workflow(args: argparse.Namespace) -> Dict[str, Any]:
    """Full cage-aware age prediction workflow."""

    # GPU ordinal control: use env var only, then pass device='cuda' to XGBoost.
    if args.device.lower() == "cuda":
        if args.gpu_id is not None and str(args.gpu_id).strip() != "":
            os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu_id)
            print(f"[GPU] CUDA_VISIBLE_DEVICES={os.environ['CUDA_VISIBLE_DEVICES']}")
        else:
            print("[GPU] device='cuda', but --gpu-id not provided. Using externally visible CUDA device(s).")
    else:
        print("[Device] CPU mode")

    print(f"[Device] XGBoost device parameter: {args.device}")

    out_path = Path(args.out)
    outdir = ensure_dir(args.outdir if args.outdir else f"{out_path.with_suffix('')}_split_tables")

    # ------------------------- Load and align -------------------------
    meta = prepare_metadata(args.meta, sample_col=args.meta_sample_col)
    ab = prepare_abundance(args.abundance, sample_col=args.abundance_sample_col)

    meta, ab = align_meta_abundance(
        meta=meta,
        ab=ab,
        target_col=args.target_col,
        group_col=args.group_col,
        mouse_id_col=args.mouse_id_col,
        diet_col=args.diet_col,
    )

    # ------------------------- Filter by diet -------------------------
    if str(args.diet_group).upper() != "ALL":
        mask = meta[args.diet_col].astype(str) == str(args.diet_group)
        meta = meta.loc[mask].copy()
        ab = ab.loc[meta.index].copy()
        if meta.empty:
            raise ValueError(f"No samples remain after filtering diet_group={args.diet_group}")

    # Clean group and target fields.
    meta[args.group_col] = meta[args.group_col].astype(str)
    meta[args.mouse_id_col] = meta[args.mouse_id_col].astype(str)
    meta[args.diet_col] = meta[args.diet_col].astype(str)
    y = meta[args.target_col].astype(float).values
    X_df = ab.copy()

    # Validate mouse nested in cage/group.
    mouse_group_counts = meta.groupby(args.mouse_id_col)[args.group_col].nunique()
    n_mice_multi_group = int((mouse_group_counts > 1).sum())
    if n_mice_multi_group > 0:
        print(
            f"[WARNING] {n_mice_multi_group} mice belong to >1 {args.group_col}. "
            "Cage-level split still blocks cage leakage, but mouse nesting is not strict."
        )

    groups = meta[args.group_col].astype(str).values
    strata_series = make_strata(
        meta=meta,
        target_col=args.target_col,
        diet_col=args.diet_col,
        age_strata_col=args.age_strata_col,
        diet_group=args.diet_group,
        n_age_bins=args.n_age_bins,
        force_diet_in_strata=args.force_diet_in_strata,
    )
    strata = strata_series.astype(str).values

    print("[Dataset]")
    print(f"  Diet group: {args.diet_group}")
    print(f"  X shape: {X_df.shape[0]:,} samples x {X_df.shape[1]:,} features")
    print(f"  Unique {args.mouse_id_col}: {meta[args.mouse_id_col].nunique():,}")
    print(f"  Unique {args.group_col}: {meta[args.group_col].nunique():,}")
    print(f"  Target: {args.target_col}")
    print(f"  Target range: {np.min(y):.3f} - {np.max(y):.3f}")
    print(f"  Stratification labels: {pd.Series(strata).nunique():,}")
    print(f"  Output directory: {outdir}")

    # ------------------------- Holdout split -------------------------
    train_idx, test_idx, holdout_report = balanced_group_holdout_split(
        y=y,
        strata=strata,
        groups=groups,
        test_size=args.test_size,
        random_state=args.random_state,
        n_candidates=args.holdout_candidates,
    )

    print("[Holdout split]")
    print(f"  Train samples: {len(train_idx):,}")
    print(f"  Test samples:  {len(test_idx):,}")
    print(f"  Train groups:  {holdout_report['n_train_groups']:,}")
    print(f"  Test groups:   {holdout_report['n_test_groups']:,}")
    print(f"  Shared groups: {holdout_report['n_shared_groups']:,}")
    print(f"  Train mean age: {holdout_report['train_mean_age']:.3f}")
    print(f"  Test mean age:  {holdout_report['test_mean_age']:.3f}")

    save_json(holdout_report, outdir / "holdout_split_report.json")
    table_paths = save_holdout_tables(
        outdir=outdir,
        X_df=X_df,
        meta=meta,
        train_idx=train_idx,
        test_idx=test_idx,
    )

    X = X_df.values.astype(np.float32)
    X_train = X[train_idx]
    X_test = X[test_idx]
    y_train = y[train_idx]
    y_test = y[test_idx]
    groups_train = groups[train_idx]
    groups_test = groups[test_idx]
    strata_train = strata[train_idx]
    strata_test = strata[test_idx]

    X_train_df = X_df.iloc[train_idx].copy()
    X_test_df = X_df.iloc[test_idx].copy()
    meta_train = meta.iloc[train_idx].copy()
    meta_test = meta.iloc[test_idx].copy()

    # ------------------------- Tuning data subset -------------------------
    if 0 < args.tune_fraction < 1.0:
        tune_train_idx_rel, _, tune_report = balanced_group_holdout_split(
            y=y_train,
            strata=strata_train,
            groups=groups_train,
            test_size=1.0 - args.tune_fraction,
            random_state=args.random_state + 101,
            n_candidates=max(500, min(args.holdout_candidates, 2000)),
        )
        X_train_tune = X_train[tune_train_idx_rel]
        y_train_tune = y_train[tune_train_idx_rel]
        groups_tune = groups_train[tune_train_idx_rel]
        strata_tune = strata_train[tune_train_idx_rel]
        save_json(tune_report, outdir / "tuning_subset_report.json")
    else:
        X_train_tune = X_train
        y_train_tune = y_train
        groups_tune = groups_train
        strata_tune = strata_train

    print("[Tuning subset]")
    print(f"  Samples: {X_train_tune.shape[0]:,}")
    print(f"  Groups:  {len(pd.unique(groups_tune)):,}")

    # ------------------------- Model params -------------------------
    base_params = {
        "subsample": args.subsample,
        "reg_lambda": args.reg_lambda,
        "reg_alpha": args.reg_alpha,
        "n_estimators": args.n_estimators,
        "min_child_weight": args.min_child_weight,
        "max_depth": args.max_depth,
        "learning_rate": args.learning_rate,
        "gamma": args.gamma,
        "colsample_bytree": args.colsample_bytree,
    }

    if args.skip_greedy:
        print("[Tuning] Skipped. Using base parameters directly.")
        used_params = dict(base_params)
        tuning_results = {
            "best_params": used_params,
            "best_params_stage1": None,
            "best_params_stage2": None,
            "stage1_cv_mae": None,
            "stage2_cv_mae": None,
            "inner_cv_name": None,
            "inner_folds": None,
            "stage1_cv_results": None,
            "stage2_cv_results": None,
        }
    else:
        tuning_results = two_stage_tuning(
            X_train_tune=X_train_tune,
            y_train_tune=y_train_tune,
            strata_tune=strata_tune,
            groups_tune=groups_tune,
            base_params=base_params,
            n_features=X_train.shape[1],
            random_state=args.random_state,
            threads=args.threads,
            device=args.device,
            inner_folds=args.inner_folds,
            n_iter_stage1=args.n_iter_stage1,
            n_iter_stage2=args.n_iter_stage2,
            search_n_jobs=args.search_n_jobs,
        )
        used_params = dict(tuning_results["best_params"])

        # Save tuning tables separately as CSV because cv_results can be large.
        if tuning_results["stage1_cv_results"] is not None:
            tuning_results["stage1_cv_results"].to_csv(outdir / "stage1_cv_results.csv", index=False)
        if tuning_results["stage2_cv_results"] is not None:
            tuning_results["stage2_cv_results"].to_csv(outdir / "stage2_cv_results.csv", index=False)

    with open(outdir / "used_xgb_params.json", "w") as f:
        json.dump(used_params, f, indent=2)

    # ------------------------- Outer CV -------------------------
    outer_splits, outer_cv_name = make_cv_splits(
        X=X_train,
        strata=strata_train,
        groups=groups_train,
        n_splits=args.outer_folds,
        random_state=args.random_state,
    )
    print(f"[Outer CV] {outer_cv_name}, folds={len(outer_splits)}")

    cv_split_records = save_cv_fold_info(
        outdir=outdir,
        X_train_df=X_train_df,
        meta_train=meta_train,
        splits=outer_splits,
        groups_train=groups_train,
        strata_train=strata_train,
        save_cv_clr_tables=args.save_cv_clr_tables,
    )

    cv_metrics: List[Dict[str, Any]] = []
    cv_pred_records: List[pd.DataFrame] = []
    cv_models: List[XGBRegressor] = [] if args.save_cv_models_in_pickle else []

    for fold_id, (tr_idx, val_idx) in enumerate(outer_splits, start=1):
        print(f"[Outer CV] Fold {fold_id}/{len(outer_splits)}")
        model = make_xgb_regressor(
            params=used_params,
            random_state=args.random_state + fold_id,
            threads=args.threads,
            device=args.device,
        )
        model.fit(X_train[tr_idx], y_train[tr_idx])
        pred = model.predict(X_train[val_idx])
        metrics = regression_metrics(y_train[val_idx], pred)
        metrics.update(
            {
                "fold": int(fold_id),
                "n_train_samples": int(len(tr_idx)),
                "n_val_samples": int(len(val_idx)),
                "n_train_groups": int(len(set(groups_train[tr_idx].astype(str)))),
                "n_val_groups": int(len(set(groups_train[val_idx].astype(str)))),
                "n_shared_groups": int(
                    len(set(groups_train[tr_idx].astype(str)).intersection(set(groups_train[val_idx].astype(str))))
                ),
            }
        )
        cv_metrics.append(metrics)
        print(f"  MAE={metrics['mae']:.4f}, RMSE={metrics['rmse']:.4f}, R2={metrics['r2']:.4f}")

        pred_df = meta_train.iloc[val_idx][
            [args.mouse_id_col, args.group_col, args.diet_col, args.target_col]
        ].copy()
        pred_df.insert(0, "sample", X_train_df.index[val_idx].astype(str))
        pred_df["fold"] = fold_id
        pred_df["strata"] = strata_train[val_idx]
        pred_df["y_true"] = y_train[val_idx]
        pred_df["y_pred"] = pred
        pred_df["abs_error"] = np.abs(y_train[val_idx] - pred)
        cv_pred_records.append(pred_df)

        if args.save_cv_models_in_pickle:
            cv_models.append(model)

    cv_metrics_df = pd.DataFrame(cv_metrics)
    cv_predictions_df = pd.concat(cv_pred_records, axis=0, ignore_index=True)
    cv_metrics_df.to_csv(outdir / "outer_cv_metrics.csv", index=False)
    cv_predictions_df.to_csv(outdir / "outer_cv_predictions.csv", index=False)

    cv_summary = {
        "outer_cv_name": outer_cv_name,
        "outer_folds": int(len(outer_splits)),
        "mae_mean": float(cv_metrics_df["mae"].mean()),
        "mae_std": float(cv_metrics_df["mae"].std(ddof=0)),
        "rmse_mean": float(cv_metrics_df["rmse"].mean()),
        "rmse_std": float(cv_metrics_df["rmse"].std(ddof=0)),
        "r2_mean": float(cv_metrics_df["r2"].mean()),
        "r2_std": float(cv_metrics_df["r2"].std(ddof=0)),
    }
    save_json(cv_summary, outdir / "outer_cv_summary.json")

    # ------------------------- Final model and holdout test -------------------------
    print("[Final model] Fitting on full holdout-training set")
    final_model = make_xgb_regressor(
        params=used_params,
        random_state=args.random_state,
        threads=args.threads,
        device=args.device,
    )
    final_model.fit(X_train, y_train)

    print("[Holdout test] Predicting held-out cage-level test set")
    y_pred_test = final_model.predict(X_test)
    holdout_metrics = regression_metrics(y_test, y_pred_test)
    print(
        f"[Holdout test] MAE={holdout_metrics['mae']:.4f}, "
        f"RMSE={holdout_metrics['rmse']:.4f}, R2={holdout_metrics['r2']:.4f}"
    )

    holdout_pred_df = meta_test[[args.mouse_id_col, args.group_col, args.diet_col, args.target_col]].copy()
    holdout_pred_df.insert(0, "sample", X_test_df.index.astype(str))
    holdout_pred_df["strata"] = strata_test
    holdout_pred_df["y_true"] = y_test
    holdout_pred_df["y_pred"] = y_pred_test
    holdout_pred_df["abs_error"] = np.abs(y_test - y_pred_test)
    holdout_pred_df.to_csv(outdir / "holdout_test_predictions.csv", index=False)
    save_json(holdout_metrics, outdir / "holdout_test_metrics.json")

    # Save model in XGBoost-native format too.
    final_model_path = outdir / "final_xgb_model.json"
    try:
        final_model.save_model(final_model_path)
    except Exception as e:
        print(f"[WARNING] Could not save native XGBoost model: {e}")
        final_model_path = None

    # ------------------------- Final results pickle -------------------------
    result: Dict[str, Any] = {
        "args": vars(args),
        "outdir": str(outdir),
        "table_paths": table_paths,
        "final_model_path": str(final_model_path) if final_model_path is not None else None,
        "final_model": final_model,
        "cv_models": cv_models if args.save_cv_models_in_pickle else None,
        "feature_names": X_df.columns.astype(str).tolist(),
        "sample_ids_all": X_df.index.astype(str).tolist(),
        "sample_ids_train": X_train_df.index.astype(str).tolist(),
        "sample_ids_test": X_test_df.index.astype(str).tolist(),
        "holdout_split_report": holdout_report,
        "cv_split_records": cv_split_records,
        "cv_metrics": cv_metrics_df,
        "cv_summary": cv_summary,
        "cv_predictions": cv_predictions_df,
        "holdout_metrics": holdout_metrics,
        "holdout_predictions": holdout_pred_df,
        "used_params_for_cv_and_final": used_params,
        "base_params": base_params,
        "best_params_stage1": tuning_results.get("best_params_stage1"),
        "best_params_stage2": tuning_results.get("best_params_stage2"),
        "stage1_cv_mae": tuning_results.get("stage1_cv_mae"),
        "stage2_cv_mae": tuning_results.get("stage2_cv_mae"),
        # Dataframes below are convenient but can make pickle large.
        "meta_train": meta_train,
        "meta_test": meta_test,
        "X_train_clr": X_train_df,
        "X_test_clr": X_test_df,
        "y_train": y_train,
        "y_test": y_test,
        "y_pred_test": y_pred_test,
    }

    # Do not keep massive cv_results DataFrames in pickle by default; saved as CSV already.
    if args.keep_tuning_cv_results_in_pickle:
        result["stage1_cv_results"] = tuning_results.get("stage1_cv_results")
        result["stage2_cv_results"] = tuning_results.get("stage2_cv_results")

    with open(out_path, "wb") as f:
        pickle.dump(result, f)

    print("[Saved]")
    print(f"  Pickle result: {out_path}")
    print(f"  Split/tables directory: {outdir}")
    print(f"  Holdout train CLR: {table_paths['holdout_train_clr']}")
    print(f"  Holdout test CLR:  {table_paths['holdout_test_clr']}")
    print(f"  CV predictions:    {outdir / 'outer_cv_predictions.csv'}")
    print(f"  Test predictions:  {outdir / 'holdout_test_predictions.csv'}")

    return result


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Cage-aware XGBoost mouse age prediction using GPU and saved split tables."
    )

    # Required I/O
    parser.add_argument("--meta", required=True, help="Metadata CSV/TSV path.")
    parser.add_argument("--abundance", required=True, help="CLR abundance CSV/TSV path.")
    parser.add_argument("--out", required=True, help="Output pickle path.")
    parser.add_argument("--outdir", default=None, help="Output directory for split tables. Default: <out>_split_tables")

    # Column names
    parser.add_argument("--meta-sample-col", default="sample", help="Sample ID column in metadata.")
    parser.add_argument("--abundance-sample-col", default=None, help="Sample ID column in abundance table. Default: infer.")
    parser.add_argument("--target-col", default="age.approx.wks", help="Regression target column.")
    parser.add_argument("--mouse-id-col", default="mouse.ID", help="Mouse ID column.")
    parser.add_argument("--group-col", default="Cage", help="Group column for leakage-safe split. Default: Cage.")
    parser.add_argument("--diet-col", default="Diet", help="Diet column.")
    parser.add_argument("--age-strata-col", default="Age bins", help="Age strata column. If missing, quantile bins are used.")

    # Analysis mode
    parser.add_argument("--diet-group", default="ALL", help="Diet group to analyze, e.g. ALL, 1D, 20, 2D, 40, AL.")
    parser.add_argument("--force-diet-in-strata", action="store_true", help="Use Diet+age strata even for diet-specific data.")
    parser.add_argument("--n-age-bins", type=int, default=10, help="Quantile age bins if --age-strata-col is unavailable.")

    # Split options
    parser.add_argument("--test-size", type=float, default=0.2, help="Held-out group-level test fraction.")
    parser.add_argument("--holdout-candidates", type=int, default=5000, help="Number of random group holdout candidates.")
    parser.add_argument("--outer-folds", type=int, default=10, help="Outer CV folds. Use 10 for ALL, usually 5 for diet-specific.")
    parser.add_argument("--inner-folds", type=int, default=5, help="Inner CV folds for tuning.")
    parser.add_argument("--random-state", type=int, default=358, help="Random seed.")

    # GPU / parallelism
    parser.add_argument("--gpu-id", default=None, help="Physical GPU ID for CUDA_VISIBLE_DEVICES, e.g. 3.")
    parser.add_argument("--device", default="cuda", choices=["cuda", "cpu"], help="XGBoost device parameter. Use cuda or cpu.")
    parser.add_argument("--threads", type=int, default=16, help="XGBoost n_jobs per model.")
    parser.add_argument("--search-n-jobs", type=int, default=1, help="RandomizedSearchCV parallel jobs. For single GPU, keep 1.")

    # Tuning control
    parser.add_argument("--skip-greedy", action="store_true", help="Skip two-stage tuning and use base params.")
    parser.add_argument("--n-iter-stage1", type=int, default=50, help="RandomizedSearchCV iterations for stage 1.")
    parser.add_argument("--n-iter-stage2", type=int, default=30, help="RandomizedSearchCV iterations for stage 2.")
    parser.add_argument("--tune-fraction", type=float, default=1.0, help="Fraction of training groups used for tuning.")

    # Base XGBoost params
    parser.add_argument("--subsample", type=float, default=0.6)
    parser.add_argument("--reg-lambda", type=float, default=2.0, dest="reg_lambda")
    parser.add_argument("--reg-alpha", type=float, default=2.0, dest="reg_alpha")
    parser.add_argument("--n-estimators", type=int, default=1200, dest="n_estimators")
    parser.add_argument("--min-child-weight", type=int, default=3, dest="min_child_weight")
    parser.add_argument("--max-depth", type=int, default=3, dest="max_depth")
    parser.add_argument("--learning-rate", type=float, default=0.01, dest="learning_rate")
    parser.add_argument("--gamma", type=float, default=1.0)
    parser.add_argument("--colsample-bytree", type=float, default=0.3, dest="colsample_bytree")

    # Save options
    parser.add_argument("--save-cv-clr-tables", action="store_true", help="Also save train/val CLR tables for every CV fold. Large output.")
    parser.add_argument("--save-cv-models-in-pickle", action="store_true", help="Store all CV models in pickle. Large output.")
    parser.add_argument("--keep-tuning-cv-results-in-pickle", action="store_true", help="Store RandomizedSearchCV result tables in pickle as well as CSV.")

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    run_workflow(args)


if __name__ == "__main__":
    main()
