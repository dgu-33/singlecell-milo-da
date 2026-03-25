# Milo Differential Abundance Analysis — CRC Single-Cell Data

Milo-based differential abundance (DA) analysis applied to colorectal cancer (CRC)
single-cell RNA-seq data, comparing MSI (microsatellite instability) vs MSS
(microsatellite stable) tumour microenvironments.

This repository contains my MILO-based differential abundance analysis module from a broader ARPA-H-related CRC single-cell analysis project.

---

## What is Milo?

[Milo](https://github.com/emdann/milopy) is a statistical method for differential
abundance testing on single-cell data. Instead of comparing predefined clusters,
it tests for abundance changes across k-nearest-neighbour (KNN) neighbourhoods,
which reduces sensitivity to arbitrary clustering decisions.

> Dann et al. (2022) *Differential abundance testing on single-cell data using
> k-nearest neighbor graphs.* Nature Biotechnology.

---

## What this script does

1. Loads a `.h5ad` single-cell dataset and maps consensus cluster labels from a CSV
2. Filters cells by microsatellite status (MSI / MSS)
3. Runs PCA, KNN graph construction, and UMAP
4. Runs the Milo DA pipeline (`make_nhoods` → `count_nhoods` → `DA_nhoods`)
5. Filters significant neighbourhoods by logFC, p-value, and spatial FDR
6. Exports MSI- and MSS-enriched cell lists as CSV files
7. Saves a DA neighbourhood graph and a volcano plot

---

## Usage

```bash
pip install -r requirements.txt

python MILO.py \
    --input_h5ad  /path/to/B_cell.h5ad \
    --input_csv   /path/to/consensus_clustering_result.csv \
    --output_dir  ./output
```

### Key arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--input_h5ad` | required | Path to the `.h5ad` input file |
| `--input_csv` | required | Path to the consensus clustering CSV |
| `--output_dir` | `./output` | Directory for output files |
| `--sample_col` | `sample_id` | obs column for biological replicates |
| `--status_col` | `microsatellite_status` | obs column for MSI/MSS condition |
| `--alpha` | `0.01` | Spatial FDR threshold |
| `--logfc_thr` | `0.2` | logFC threshold |
| `--pval_thr` | `0.05` | p-value threshold |
| `--seed` | `42` | Random seed |

---

## Project structure

```
milo-crc-da/
├── MILO.py          # DA analysis pipeline
├── requirements.txt
├── .gitignore
└── README.md
```

---

## Note on data

Input data files (`.h5ad`, `.csv`) are not included in this repository as they
are part of a larger internal project dataset and are too large for GitHub.
