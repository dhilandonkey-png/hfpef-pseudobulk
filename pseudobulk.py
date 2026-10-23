"""Turn single-nucleus counts into pseudobulk samples (one per donor x cell type).

Follows Hahn, Chaffin et al., Circ Res 2026 (HFpEF snRNA-seq, SCP3342):
  * sum raw UMI counts of all nuclei of a cell type within each patient
  * keep a cell type only if >= 20 nuclei in >= 3 individuals from BOTH groups
  * also build an "all cell types" pseudobulk (one sample per patient)

Example:
  python pseudobulk.py --h5ad SCP3342.h5ad --out results/ \\
      --donor-col donor_id --celltype-col cell_type --group-col disease \\
      --covariates sex pool
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from load_data import _looks_like_counts, load


def aggregate(adata: ad.AnnData, keys: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Sum counts over cells sharing the same values of `keys`.

    Returns (counts: samples x genes, sample metadata with n_cells).
    """
    groups = adata.obs[keys].astype(str).agg("|".join, axis=1)
    codes, labels = pd.factorize(groups, sort=True)
    # one-hot indicator (samples x cells) @ counts (cells x genes) = summed counts
    indicator = sp.csr_matrix(
        (np.ones(len(codes)), (codes, np.arange(len(codes)))),
        shape=(len(labels), len(codes)),
    )
    summed = indicator @ adata.X
    summed = summed.toarray() if sp.issparse(summed) else np.asarray(summed)
    counts = pd.DataFrame(np.rint(summed).astype(np.int64), index=labels, columns=adata.var_names)

    meta = adata.obs[keys].astype(str).copy()
    meta.index = groups.values
    meta = meta[~meta.index.duplicated()].loc[labels]
    meta["n_cells"] = np.bincount(codes)
    return counts, meta


def eligible_celltypes(
    meta: pd.DataFrame, celltype_col: str, group_col: str, min_cells: int, min_donors: int
) -> list[str]:
    """Cell types with >= min_cells nuclei in >= min_donors donors of every group."""
    keep = []
    n_groups = meta[group_col].nunique()
    for ct, sub in meta.groupby(celltype_col):
        ok = sub[sub["n_cells"] >= min_cells].groupby(group_col).size()
        if len(ok) == n_groups and (ok >= min_donors).all():
            keep.append(ct)
    return keep


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = p.add_argument_group("input")
    src.add_argument("--h5ad")
    src.add_argument("--layer", help="AnnData layer holding raw counts (e.g. 'counts')")
    src.add_argument("--matrix")
    src.add_argument("--barcodes")
    src.add_argument("--features")
    src.add_argument("--metadata", help="SCP metadata.txt (NAME + TYPE rows)")
    p.add_argument("--out", required=True)
    p.add_argument("--donor-col", required=True, help="patient / donor id column")
    p.add_argument("--celltype-col", required=True)
    p.add_argument("--group-col", required=True, help="condition column, e.g. HFpEF vs control")
    p.add_argument("--covariates", nargs="*", default=[], help="donor-level columns kept in metadata (sex, pool, ...)")
    p.add_argument("--min-cells", type=int, default=20, help="min nuclei per donor x cell type (paper: 20)")
    p.add_argument("--min-donors", type=int, default=3, help="min eligible donors per group (paper: 3)")
    p.add_argument("--allow-non-integer", action="store_true", help="aggregate even if X is not raw counts")
    args = p.parse_args(argv)

    adata = load(args.h5ad, args.matrix, args.barcodes, args.features, args.metadata, args.layer)
    print(f"loaded {adata.n_obs:,} nuclei x {adata.n_vars:,} genes")

    needed = [args.donor_col, args.celltype_col, args.group_col, *args.covariates]
    missing = [c for c in needed if c not in adata.obs.columns]
    if missing:
        raise SystemExit(f"missing obs columns {missing}; available: {list(adata.obs.columns)}")
    if not _looks_like_counts(adata.X) and not args.allow_non_integer:
        raise SystemExit(
            "X does not look like raw integer counts. Pseudobulk must sum RAW counts; "
            "pass --layer counts, or --allow-non-integer if you know what you are doing."
        )
    adata = adata[adata.obs[needed].notna().all(axis=1)].copy()

    # donor-level covariates must be constant within a donor
    donor_meta = adata.obs[[args.donor_col, args.group_col, *args.covariates]].astype(str).drop_duplicates()
    if donor_meta[args.donor_col].duplicated().any():
        bad = donor_meta[donor_meta[args.donor_col].duplicated(keep=False)]
        raise SystemExit(
            f"group/covariates vary within a donor:\n{bad}\n"
            f"Is --donor-col '{args.donor_col}' really one patient? Pool/library columns "
            "(e.g. SCP3342's biosample_id = pool2..pool10) mix HFpEF and control patients."
        )

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    # 1) per donor x cell type
    counts, meta = aggregate(adata, [args.donor_col, args.celltype_col])
    meta = meta.merge(donor_meta, on=args.donor_col, how="left").set_index(counts.index)
    keep = eligible_celltypes(meta, args.celltype_col, args.group_col, args.min_cells, args.min_donors)
    print(f"cell types eligible for DE ({len(keep)}): {keep}")

    summary = (
        meta.pivot_table(index=args.celltype_col, columns=args.group_col, values="n_cells", aggfunc="sum")
        .fillna(0)
        .astype(int)
    )
    summary["eligible_for_DE"] = summary.index.isin(keep)
    summary.to_csv(out / "celltype_summary.csv")

    ct_dir = out / "by_celltype"
    ct_dir.mkdir(exist_ok=True)
    for ct in keep:
        rows = (meta[args.celltype_col] == ct) & (meta["n_cells"] >= args.min_cells)
        safe = "".join(ch if ch.isalnum() else "_" for ch in ct)
        sub_meta = meta.loc[rows].set_index(args.donor_col)
        sub_counts = counts.loc[rows].set_axis(sub_meta.index)
        # genes x samples, the orientation limma / edgeR / DESeq2 expect
        sub_counts.T.to_csv(ct_dir / f"{safe}_counts.csv.gz")
        sub_meta.to_csv(ct_dir / f"{safe}_metadata.csv")

    ad.AnnData(X=sp.csr_matrix(counts.values), obs=meta, var=pd.DataFrame(index=counts.columns)).write_h5ad(
        out / "pseudobulk_donor_x_celltype.h5ad"
    )

    # 2) all cell types combined, one sample per donor
    all_counts, all_meta = aggregate(adata, [args.donor_col])
    all_meta = all_meta.merge(donor_meta, on=args.donor_col, how="left").set_index(args.donor_col)
    all_counts.index = all_meta.index
    all_counts.T.to_csv(out / "all_celltypes_counts.csv.gz")
    all_meta.to_csv(out / "all_celltypes_metadata.csv")

    print(f"wrote {len(counts)} donor x cell-type pseudobulk samples to {out}/")


if __name__ == "__main__":
    main()
