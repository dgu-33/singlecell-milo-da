#!/usr/bin/env python3
"""
MILO.py — Milo Differential Abundance Analysis for CRC Single-Cell Data

Applies Milo-based differential abundance (DA) analysis to single-cell RNA-seq
data from colorectal cancer (CRC) samples, comparing MSI (microsatellite
instable) vs MSS (microsatellite stable) tumour microenvironments.

This script covers the DA analysis step of a larger CRC atlas pipeline.
Upstream steps (cell type annotation, consensus clustering) were produced by
the broader team; this script takes their outputs as inputs.

Usage:
    python MILO.py \\
        --input_h5ad  /path/to/B_cell.h5ad \\
        --input_csv   /path/to/consensus_clustering_result.csv \\
        --output_dir  ./output

Reference:
    Dann et al. (2022) Differential abundance testing on single-cell data
    using k-nearest neighbor graphs. Nature Biotechnology.
    https://github.com/emdann/milopy
"""

import os
import argparse

import h5py
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix

from milopy.core import make_nhoods, count_nhoods, DA_nhoods
from milopy.utils import build_nhood_graph
from milopy.plot import plot_nhood_graph


# =============================================================================
# 1. Data loading
# =============================================================================

def load_h5ad(file_path: str) -> anndata.AnnData:
    """Load an AnnData object from an .h5ad file using h5py.

    A custom h5py-based loader is used here to handle non-standard or
    partially-written .h5ad files that anndata.read_h5ad() may reject.
    Handles both dense and sparse X matrices, categorical obs fields,
    and byte-string decoding.

    Args:
        file_path: Full path to the .h5ad file.

    Returns:
        AnnData with total-count normalised and log1p-transformed expression.
    """
    with h5py.File(file_path, 'r', libver='latest') as f:
        X_sparse = _read_X(f)
        obs_df = _read_obs(f)
        var_names = _read_var_names(f)

    adata = anndata.AnnData(X=X_sparse, obs=obs_df)
    adata.var_names = var_names
    adata.var = pd.DataFrame(index=adata.var_names)
    adata.raw = adata.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def _read_X(f: h5py.File) -> csr_matrix:
    """Parse the X matrix from an h5py file handle."""
    X_obj = f['X']
    if isinstance(X_obj, h5py.Dataset):
        return csr_matrix(X_obj[()])
    elif isinstance(X_obj, h5py.Group):
        shape = X_obj.attrs.get('shape')
        if shape is None:
            var_names = _read_var_names(f)
            n_cells = len(X_obj['indptr'][()]) - 1
            shape = (n_cells, len(var_names))
        return csr_matrix(
            (X_obj['data'][()], X_obj['indices'][()], X_obj['indptr'][()]),
            shape=tuple(shape)
        )
    raise TypeError(f"Unexpected type for X matrix: {type(X_obj)}")


def _read_obs(f: h5py.File) -> pd.DataFrame:
    """Parse the obs DataFrame from an h5py file handle."""
    obs_index_ds = f.get('obs/_index')
    if obs_index_ds is None:
        raise KeyError("'obs/_index' not found in file.")
    obs_index = [s.decode('utf-8') if isinstance(s, bytes) else s for s in obs_index_ds[()]]

    obs_group = f.get('obs')
    obs_data = {}
    if obs_group is not None:
        for key in obs_group.keys():
            item = obs_group.get(key)
            if isinstance(item, h5py.Group) and "codes" in item and "categories" in item:
                cats = [c.decode("utf-8") if isinstance(c, bytes) else c
                        for c in item["categories"][()]]
                obs_data[key] = [cats[c] for c in item["codes"][()]]
            elif isinstance(item, h5py.Group) and "codes" in item:
                obs_data[key] = [c.decode('utf-8') if isinstance(c, bytes) else c
                                 for c in item["codes"][()]]
            elif isinstance(item, h5py.Dataset):
                val = item[()]
                if val.dtype.type is np.bytes_:
                    val = [v.decode('utf-8') for v in val]
                obs_data[key] = val

    return pd.DataFrame(obs_data, index=obs_index)


def _read_var_names(f: h5py.File) -> list:
    """Parse variable (gene) names from an h5py file handle."""
    var_group = f.get('var')
    if var_group is None:
        raise KeyError("'var' group not found in file.")
    ds = var_group.get('var_names') or var_group.get('_index')
    if ds is None:
        raise KeyError("Neither 'var/var_names' nor 'var/_index' found.")
    return [s.decode('utf-8') if isinstance(s, bytes) else s for s in ds[()]]


# =============================================================================
# 2. Preprocessing and annotation mapping
# =============================================================================

def map_consensus_clusters(adata: anndata.AnnData, clust_df: pd.DataFrame,
                            cell_type_col: str = 'cell_type') -> anndata.AnnData:
    """Map consensus cluster labels from the clustering CSV onto adata.obs.

    The clustering CSV is indexed by cell name, with cell-type names as columns
    and cluster IDs as values.

    Args:
        adata: AnnData object whose obs contains cell_type_col.
        clust_df: DataFrame indexed by cell name with cell-type columns.
        cell_type_col: obs column holding cell type labels.

    Returns:
        AnnData with 'consensus_cluster' added to obs, filtered to cells
        that have a valid cluster assignment.
    """
    def _lookup(row):
        cell_type = row.get(cell_type_col)
        if (row.name in clust_df.index
                and pd.notnull(cell_type)
                and cell_type in clust_df.columns):
            return clust_df.loc[row.name, cell_type]
        return np.nan

    adata.obs['consensus_cluster'] = adata.obs.apply(_lookup, axis=1)
    adata = adata[adata.obs['consensus_cluster'].notna()].copy()
    print(f"Cells after consensus cluster filter: {adata.n_obs}")
    return adata


def filter_by_msi_status(adata: anndata.AnnData,
                          status_col: str = 'microsatellite_status') -> anndata.AnnData:
    """Normalise and filter cells by microsatellite status (MSI / MSS).

    Maps 'MSI-H' to 'MSI' and removes cells with unknown/missing status.

    Args:
        adata: AnnData object.
        status_col: obs column containing microsatellite status.

    Returns:
        Filtered AnnData with status as an ordered categorical (MSS < MSI).
    """
    mapping = {"MSI": "MSI", "MSI-H": "MSI", "MSS": "MSS"}
    adata.obs[status_col] = adata.obs[status_col].astype(str).map(mapping)
    adata = adata[adata.obs[status_col].notna()].copy()
    adata.obs[status_col] = pd.Categorical(
        adata.obs[status_col], categories=["MSS", "MSI"], ordered=True
    )
    print(f"Microsatellite status distribution:\n{adata.obs[status_col].value_counts().to_string()}")
    return adata


def run_preprocessing(adata: anndata.AnnData,
                       n_pcs: int = 30,
                       n_neighbors: int = 15) -> anndata.AnnData:
    """Run PCA, KNN graph construction, and UMAP embedding.

    Args:
        adata: Normalised AnnData.
        n_pcs: Number of principal components.
        n_neighbors: Number of neighbors for the KNN graph.

    Returns:
        AnnData with X_pca, neighbors, and X_umap populated.
    """
    sc.pp.pca(adata, n_comps=n_pcs)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata)
    return adata


# =============================================================================
# 3. Milo DA analysis
# =============================================================================

def run_milo(adata: anndata.AnnData, sample_col: str, design: str,
             prop: float = 0.1, seed: int = 42) -> anndata.AnnData:
    """Run the full Milo differential abundance pipeline.

    Steps: make_nhoods -> count_nhoods -> DA_nhoods -> build_nhood_graph.

    Args:
        adata: Preprocessed AnnData with UMAP coordinates.
        sample_col: obs column identifying biological replicates/samples.
        design: R-style formula for the DA test (e.g. '~ microsatellite_status').
        prop: Proportion of cells to use as neighbourhood index points.
        seed: Random seed for reproducibility.

    Returns:
        AnnData with Milo results in adata.uns['nhood_adata'].
    """
    make_nhoods(adata, neighbors_key=None, prop=prop, seed=seed)
    count_nhoods(adata, sample_col=sample_col)
    DA_nhoods(adata, design=design)
    build_nhood_graph(adata, basis="X_umap")
    return adata


# =============================================================================
# 4. Result filtering and cell export
# =============================================================================

def filter_neighborhoods(df: pd.DataFrame, logfc_thr: float, pval_thr: float,
                          fdr_thr: float, sign: str) -> pd.DataFrame:
    """Filter Milo neighbourhood results by significance thresholds.

    Args:
        df: nhood_adata.obs with 'logFC', 'PValue', 'SpatialFDR' columns.
        logfc_thr: Minimum absolute logFC.
        pval_thr: Maximum raw p-value.
        fdr_thr: Maximum spatial FDR.
        sign: 'pos' for MSI-enriched neighbourhoods, 'neg' for MSS-enriched.

    Returns:
        Filtered DataFrame of significant neighbourhoods.
    """
    cond = (df["PValue"] < pval_thr) & (df["SpatialFDR"] < fdr_thr)
    if sign == "pos":
        cond = cond & (df["logFC"] >= logfc_thr)
    elif sign == "neg":
        cond = cond & (df["logFC"] <= -logfc_thr)
    else:
        raise ValueError("sign must be 'pos' or 'neg'.")
    return df[cond]


def get_cells_in_nhoods(nhood_names_list, nhoods_mat,
                         adata_obs_names, nhood_index: list) -> set:
    """Collect all cell names belonging to a set of neighbourhoods.

    Args:
        nhood_names_list: Neighbourhood names to collect cells from.
        nhoods_mat: Cell x neighbourhood membership matrix (adata.obsm['nhoods']).
        adata_obs_names: Cell names (adata.obs_names).
        nhood_index: Ordered list of all neighbourhood names.

    Returns:
        Set of cell names.
    """
    cells = set()
    for nh_name in nhood_names_list:
        col_idx = nhood_index.index(nh_name)
        membership = nhoods_mat[:, col_idx]
        if hasattr(membership, "toarray"):
            membership = membership.toarray().ravel()
        cells.update(adata_obs_names[np.where(membership > 0)[0]])
    return cells


def export_cell_lists(msi_cells: set, mss_cells: set, output_dir: str) -> None:
    """Save MSI- and MSS-enriched cell lists to CSV files.

    Args:
        msi_cells: Cell names enriched in MSI neighbourhoods.
        mss_cells: Cell names enriched in MSS neighbourhoods.
        output_dir: Directory for output CSV files.
    """
    pd.DataFrame({"cell_name": sorted(msi_cells)}).to_csv(
        os.path.join(output_dir, "msi_enriched_cells.csv"), index=False
    )
    pd.DataFrame({"cell_name": sorted(mss_cells)}).to_csv(
        os.path.join(output_dir, "mss_enriched_cells.csv"), index=False
    )
    print(f"MSI-enriched cells ({len(msi_cells)}) -> msi_enriched_cells.csv")
    print(f"MSS-enriched cells ({len(mss_cells)}) -> mss_enriched_cells.csv")


# =============================================================================
# 5. Visualisation
# =============================================================================

def plot_da_graph(adata: anndata.AnnData, output_dir: str,
                  alpha: float = 0.01) -> None:
    """Plot the Milo DA neighbourhood graph coloured by logFC and save to disk."""
    plot_nhood_graph(
        adata,
        alpha=alpha,
        min_logFC=0,
        min_size=2,
        plot_edges=False,
        title="DA log-Fold Change"
    )
    out_path = os.path.join(output_dir, "DA_logFC.png")
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"DA graph saved: {out_path}")


def plot_volcano(nhood_df: pd.DataFrame, output_dir: str,
                 logfc_thr: float = 0.2, pval_thr: float = 0.05,
                 fdr_thr: float = 0.01) -> None:
    """Plot a volcano plot of logFC vs -log10(p-value) for all neighbourhoods.

    Colours: red = MSI-enriched, blue = MSS-enriched, gray = not significant.

    Args:
        nhood_df: nhood_adata.obs DataFrame containing Milo results.
        output_dir: Directory to save the plot.
        logfc_thr: logFC threshold lines.
        pval_thr: p-value threshold line.
        fdr_thr: FDR threshold used for significance colouring.
    """
    df = nhood_df.copy()
    df["minus_log10_pval"] = -np.log10(df["PValue"] + 1e-10)

    def _status(row):
        if row["PValue"] < pval_thr and row["SpatialFDR"] < fdr_thr:
            if row["logFC"] >= logfc_thr:
                return "MSI high"
            elif row["logFC"] <= -logfc_thr:
                return "MSS high"
        return "Not significant"

    df["status"] = df.apply(_status, axis=1)
    color_map = {"MSI high": "red", "MSS high": "blue", "Not significant": "gray"}

    fig, ax = plt.subplots(figsize=(8, 6))
    for status, group in df.groupby("status"):
        ax.scatter(group["logFC"], group["minus_log10_pval"],
                   s=10, c=color_map[status], alpha=0.7, label=status)

    ax.axvline(x=logfc_thr,  color="red",   linestyle="--", linewidth=1)
    ax.axvline(x=-logfc_thr, color="blue",  linestyle="--", linewidth=1)
    ax.axhline(y=-np.log10(pval_thr), color="green", linestyle="--", linewidth=1,
               label=f"p < {pval_thr}")
    ax.set_xlabel("log Fold Change")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title("Volcano Plot — Milo Neighbourhood DA (MSI vs MSS)")
    ax.legend()
    plt.tight_layout()

    out_path = os.path.join(output_dir, "volcano_plot.png")
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Volcano plot saved: {out_path}")


# =============================================================================
# 6. CLI
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Milo differential abundance analysis on CRC single-cell data."
    )
    p.add_argument("--input_h5ad",     required=True,
                   help="Path to the input .h5ad file.")
    p.add_argument("--input_csv",      required=True,
                   help="Path to the consensus clustering CSV (indexed by cell_name).")
    p.add_argument("--output_dir",     default="./output",
                   help="Directory for output files. Created if it does not exist.")
    p.add_argument("--sample_col",     default="sample_id",
                   help="obs column identifying biological samples/replicates.")
    p.add_argument("--status_col",     default="microsatellite_status",
                   help="obs column for condition labels (MSI/MSS).")
    p.add_argument("--cell_type_col",  default="cell_type",
                   help="obs column for cell type labels.")
    p.add_argument("--design",         default="~ microsatellite_status",
                   help="Design formula for the Milo DA test.")
    p.add_argument("--prop",           type=float, default=0.1,
                   help="Proportion of cells used as neighbourhood index points.")
    p.add_argument("--alpha",          type=float, default=0.01,
                   help="FDR threshold for significance.")
    p.add_argument("--logfc_thr",      type=float, default=0.2,
                   help="logFC threshold for filtering significant neighbourhoods.")
    p.add_argument("--pval_thr",       type=float, default=0.05,
                   help="p-value threshold for filtering.")
    p.add_argument("--seed",           type=int,   default=42,
                   help="Random seed for reproducibility.")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # --- Load data ---
    print(f"Loading: {args.input_h5ad}")
    adata = load_h5ad(args.input_h5ad)

    print(f"Loading clustering CSV: {args.input_csv}")
    clust_df = pd.read_csv(args.input_csv).set_index('cell_name')

    # --- Annotation mapping and filtering ---
    adata = map_consensus_clusters(adata, clust_df, cell_type_col=args.cell_type_col)
    adata = filter_by_msi_status(adata, status_col=args.status_col)

    # --- Preprocessing ---
    print("Running PCA, neighbors, UMAP...")
    adata = run_preprocessing(adata)

    # --- Milo DA analysis ---
    print("Running Milo DA analysis...")
    adata = run_milo(adata, sample_col=args.sample_col, design=args.design,
                     prop=args.prop, seed=args.seed)

    nhood_df   = adata.uns["nhood_adata"].obs.copy()
    nhood_names = nhood_df.index.tolist()
    nhoods_mat  = adata.obsm["nhoods"]

    print(f"Total neighbourhoods: {len(nhood_df)}")
    print(f"Significant (SpatialFDR < {args.alpha}): "
          f"{(nhood_df['SpatialFDR'] < args.alpha).sum()}")

    # --- Filter and export significant neighbourhoods ---
    msi_nhoods = filter_neighborhoods(nhood_df, args.logfc_thr, args.pval_thr, args.alpha, "pos")
    mss_nhoods = filter_neighborhoods(nhood_df, args.logfc_thr, args.pval_thr, args.alpha, "neg")
    print(f"MSI-enriched neighbourhoods: {len(msi_nhoods)}")
    print(f"MSS-enriched neighbourhoods: {len(mss_nhoods)}")

    msi_cells = get_cells_in_nhoods(msi_nhoods.index, nhoods_mat, adata.obs_names, nhood_names)
    mss_cells = get_cells_in_nhoods(mss_nhoods.index, nhoods_mat, adata.obs_names, nhood_names)
    export_cell_lists(msi_cells, mss_cells, args.output_dir)

    # --- Plots ---
    plot_da_graph(adata, args.output_dir, alpha=args.alpha)
    plot_volcano(nhood_df, args.output_dir,
                 logfc_thr=args.logfc_thr, pval_thr=args.pval_thr, fdr_thr=args.alpha)

    print("Done.")


if __name__ == "__main__":
    main()
