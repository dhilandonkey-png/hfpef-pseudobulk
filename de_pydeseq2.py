"""Python-only pseudobulk DE with PyDESeq2 (alternative to limma_voom.R).

The paper used limma-voom; DESeq2 is a different model, so gene lists will
overlap but not match exactly. Use this when R is not available.

  python de_pydeseq2.py results/ --group-col disease --case HFpEF --control control \\
      --covariates sex pool
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats


def run_de(counts_file: Path, meta_file: Path, group_col: str, case: str, control: str, covariates: list[str]):
    counts = pd.read_csv(counts_file, index_col=0).T  # samples x genes
    meta = pd.read_csv(meta_file, index_col=0).loc[counts.index]
    covs = [c for c in covariates if meta[c].nunique() > 1]
    meta = meta[[group_col, *covs]].astype(str)
    meta = meta.rename(columns={group_col: "group"})

    # simple expression filter: >= 10 counts in at least the size of the smaller group
    min_n = meta["group"].value_counts().min()
    counts = counts.loc[:, (counts >= 10).sum(axis=0) >= min_n]

    dds = DeseqDataSet(
        counts=counts,
        metadata=meta,
        design="~" + " + ".join([*covs, "group"]),
        refit_cooks=True,
        inference=DefaultInference(n_cpus=1),
        quiet=True,
    )
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["group", case, control], inference=DefaultInference(n_cpus=1), quiet=True)
    stats.summary()
    return stats.results_df.sort_values("pvalue")


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("pseudobulk_dir")
    p.add_argument("--group-col", required=True)
    p.add_argument("--case", required=True)
    p.add_argument("--control", required=True)
    p.add_argument("--covariates", nargs="*", default=[])
    args = p.parse_args(argv)

    root = Path(args.pseudobulk_dir)
    out = root / "de_pydeseq2"
    out.mkdir(exist_ok=True)
    jobs = [
        (f.name.removesuffix("_counts.csv.gz"), f, f.with_name(f.name.replace("_counts.csv.gz", "_metadata.csv")))
        for f in sorted((root / "by_celltype").glob("*_counts.csv.gz"))
    ]
    jobs.append(("all_celltypes", root / "all_celltypes_counts.csv.gz", root / "all_celltypes_metadata.csv"))

    for label, counts_file, meta_file in jobs:
        res = run_de(counts_file, meta_file, args.group_col, args.case, args.control, args.covariates)
        res.to_csv(out / f"{label}_de.csv")
        print(f"{label:25s} genes={len(res):6d} DE(padj<0.05)={(res['padj'] < 0.05).sum()}")


if __name__ == "__main__":
    main()
