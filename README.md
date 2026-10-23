# HFpEF snRNA-seq → pseudobulk

Pipeline for turning the single-nucleus data from
**Hahn, Chaffin et al., "Single Cell Analysis of Human Heart Failure with Preserved Ejection Fraction", Circ Res 2026**
([PMC13286016](https://pmc.ncbi.nlm.nih.gov/articles/PMC13286016/), doi:10.1161/CIRCRESAHA.125.327433)
into pseudobulk samples, then running differential expression the way the paper did.

## The dataset

| | |
|---|---|
| Assay | **single-nucleus** RNA-seq (snRNA-seq), 10x Chromium 3' v3.1 |
| Tissue | human septal endomyocardial biopsies |
| Donors | 19 HFpEF, 24 non-failing controls (pooled 6 per library, demultiplexed with souporcell) |
| Nuclei after QC | 48,866, 14 cell types |
| Upstream | CellRanger v4 (pre-mRNA GRCh38-2020-A) → CellBender → souporcell → scanpy QC |
| Data | Broad Single Cell Portal **SCP3342**: https://singlecell.broadinstitute.org/single_cell/study/SCP3342 |

## What "pseudobulk" means here

snRNA-seq gives thousands of nuclei per patient. Nuclei from the same person are not
independent, so testing each nucleus as a replicate inflates false positives
("pseudoreplication"). Pseudobulk fixes that:

```
nuclei x genes (raw UMI counts)
        │  group by (patient, cell type), SUM counts
        ▼
(patient x cell type) x genes   ← each row looks like a small bulk RNA-seq sample
        │  per cell type: filter genes, TMM normalise, voom, limma ~ disease + sex + pool
        ▼
DE table per cell type (logFC, P, BH-adjusted P)
```

Rules taken from the paper:
- sum **raw counts** (not normalised / log values)
- a cell type is tested only if it has **≥ 20 nuclei in ≥ 3 individuals from both groups**
- model adjusted for **sex** and **pool**; Benjamini-Hochberg FDR 0.05 within each cell type
- an extra "all cell types" pseudobulk (one sample per patient) for comparison with bulk RNA-seq

## 1. Download the data

1. Sign in at the SCP3342 link above (the portal requires a Google login).
2. Open **Download** and get the **"snRNA-seq anndata file, h5ad"** (48,866 nuclei x 36,601 genes).
3. List its columns and check it holds raw counts:
   ```bash
   python -I inspect_h5ad.py SCP3342.h5ad
   ```
   Pick the columns for **patient**, **cell type**, **disease group** and **sex**.

> **`biosample_id` is the pool, not the patient.** On SCP3342 it takes 8 values
> (pool2, 3, 4, 5, 7, 8, 9, 10; pools 1, 6 and 6b failed QC in the paper). Each pool mixes
> 3 HFpEF and 3 control patients, so use it as the `pool` **covariate**, never as `--donor-col`.
> The patient column comes from souporcell demultiplexing and should have 43 values.

Make sure the matrix holds **raw counts**. If the portal only provides normalised values, the script stops and tells you so.

## Quick start for the SCP3342 file

The portal's `HFpEF_snRNAseq_single_cell_portal_withRawCounts_09.24.2026.h5ad` (48,866 nuclei x 36,601 genes) has:

| role | column / layer |
|---|---|
| patient | `donor_id` (43) |
| cell type | `original_cell_type` (the paper's 14 labels) |
| group | `original_disease` (`hfpef` / `control`) |
| covariates | `sex`, `biosample_id` (= pool, 8 values) |
| raw counts | layers `raw_count_cellbender` (main analysis) and `raw_count_cellranger` (paper's concordance check) |

`X` is normalised, so always pass `--layer`.

```bash
H5AD=~/Downloads/HFpEF_snRNAseq_single_cell_portal_withRawCounts_09.24.2026.h5ad

python pseudobulk.py --h5ad "$H5AD" --layer raw_count_cellbender --out results_cellbender \
    --donor-col donor_id --celltype-col original_cell_type --group-col original_disease \
    --covariates sex biosample_id

python de_pydeseq2.py results_cellbender --group-col original_disease \
    --case hfpef --control control --covariates sex biosample_id
# or, as in the paper:
Rscript limma_voom.R results_cellbender original_disease hfpef control sex biosample_id
```

Repeat with `--layer raw_count_cellranger --out results_cellranger` to reproduce the paper's
requirement that a gene be significant in the same direction in both matrices.

## 2. Install

```bash
git clone https://github.com/dhilandonkey-png/hfpef-pseudobulk.git
cd hfpef-pseudobulk
pip install -r requirements.txt
# for the paper's exact DE method (R):
#   install.packages("BiocManager"); BiocManager::install(c("limma", "edgeR"))
```

## 3. Build pseudobulk

```bash
# from an h5ad (use --layer counts if raw counts live in a layer)
python pseudobulk.py --h5ad SCP3342.h5ad --out results \
    --donor-col <donor> --celltype-col <cell_type> --group-col <disease> \
    --covariates <sex> <pool>

# or from SCP's matrix + metadata files
python pseudobulk.py --matrix matrix.mtx.gz --barcodes barcodes.tsv.gz \
    --features features.tsv.gz --metadata metadata.txt --out results \
    --donor-col <donor> --celltype-col <cell_type> --group-col <disease> \
    --covariates <sex> <pool>
```

Outputs in `results/`:

| file | contents |
|---|---|
| `celltype_summary.csv` | nuclei per cell type per group, and whether it passes the DE filter |
| `by_celltype/<CT>_counts.csv.gz` | genes x patients summed counts |
| `by_celltype/<CT>_metadata.csv` | patient, group, sex, pool, n_cells |
| `all_celltypes_counts.csv.gz` / `_metadata.csv` | one pseudobulk sample per patient |
| `pseudobulk_donor_x_celltype.h5ad` | everything as one AnnData |

## 4. Differential expression

Paper method (limma-voom, R):
```bash
Rscript limma_voom.R results <disease_col> HFpEF <control_label> <sex_col> <pool_col>
```

Python alternative (PyDESeq2, a different model, so results will be close but not identical):
```bash
python de_pydeseq2.py results --group-col <disease_col> --case HFpEF --control <control_label> \
    --covariates <sex_col> <pool_col>
```

## Try it without the real data

```bash
python make_synthetic.py synth.h5ad   # fake data, 43 donors, 100 planted DE genes in cardiomyocytes
python pseudobulk.py --h5ad synth.h5ad --out synth_res --donor-col donor_id \
    --celltype-col cell_type --group-col disease --covariates sex pool
python de_pydeseq2.py synth_res --group-col disease --case HFpEF --control control --covariates sex pool
```

On that synthetic set the pipeline drops the rare `Lymphocyte` type (fails the
≥20-nuclei rule). In cardiomyocytes PyDESeq2 recovers all 77 planted genes that
pass the expression filter, plus 8 false positives. The other cell types, which
have no planted signal, return 0 DE genes.

## Not reproduced here

- Upstream processing (CellRanger, CellBender, souporcell, Scrublet, QC). The SCP files are already processed.
- The paper's extra significance rule: concordant in **both** CellBender and CellRanger matrices, plus a background-contamination heuristic ≤ 0.40. This needs both matrices and the supplement's definition.
- Exact gene expression thresholds, which are in the paper's Data Supplement PDF. `limma_voom.R` uses edgeR's `filterByExpr` defaults.
