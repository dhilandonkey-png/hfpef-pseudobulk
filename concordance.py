"""Keep genes significant in the same direction with CellBender AND CellRanger counts.

This is the paper's significance rule (Hahn, Chaffin et al., Circ Res 2026):
adjusted P < 0.05 in both matrices, with the same sign of log2 fold change.
(The paper's extra background-contamination heuristic <= 0.40 is not applied.)
Also draws a volcano plot per cell type in the style of the paper's Figure 2a.

  python concordance.py results_cellbender/de_pydeseq2 results_cellranger/de_pydeseq2 --out final
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

UP, DOWN, NS = "#e34948", "#2a78d6", "#c3c2b7"  # diverging red / blue, gray = not significant
TEXT, MUTED, SURFACE = "#0b0b0b", "#52514e", "#fcfcfb"


def concordant(cb: pd.DataFrame, cr: pd.DataFrame, alpha: float) -> pd.DataFrame:
    cb = cb[["baseMean", "log2FoldChange", "pvalue", "padj"]].rename(columns={"log2FoldChange": "log2FC"})
    cr = cr[["log2FoldChange", "padj"]].rename(columns={"log2FoldChange": "log2FC_cellranger", "padj": "padj_cellranger"})
    df = cb.join(cr, how="inner")
    df["significant"] = (
        (df["padj"] < alpha)
        & (df["padj_cellranger"] < alpha)
        & (np.sign(df["log2FC"]) == np.sign(df["log2FC_cellranger"]))
    )
    df["direction"] = np.where(df["significant"], np.where(df["log2FC"] > 0, "up", "down"), "")
    return df.sort_values("pvalue")


def volcano(df: pd.DataFrame, title: str, path: Path, n_labels: int = 10) -> None:
    df = df.dropna(subset=["log2FC", "pvalue"])
    y = -np.log10(df["pvalue"].clip(lower=1e-300))
    fig, ax = plt.subplots(figsize=(6, 5), dpi=150, facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    for direction, color, label in [("", NS, "not significant"), ("down", DOWN, "down in HFpEF"), ("up", UP, "up in HFpEF")]:
        m = df["direction"] == direction
        ax.scatter(df.loc[m, "log2FC"], y[m], s=8 if direction else 4, c=color, linewidths=0,
                   label=f"{label} ({m.sum():,})", rasterized=True)
    # label the top genes each way, skipping any that would overlap a label already placed
    ax.autoscale_view()
    placed: list[np.ndarray] = []
    for direction in ("up", "down"):
        top = df[df["direction"] == direction].head(n_labels)
        for gene, row in top.iterrows():
            xy = (row["log2FC"], -np.log10(max(row["pvalue"], 1e-300)))
            px = ax.transData.transform(xy)
            if any(abs(px[0] - q[0]) < 60 and abs(px[1] - q[1]) < 14 for q in placed):
                continue
            placed.append(px)
            ax.annotate(gene, xy, fontsize=6, color=TEXT, xytext=(3, 2), textcoords="offset points")
    ax.axvline(0, color=MUTED, lw=0.5)
    ax.set_xlabel("log2 fold change (HFpEF vs control, CellBender)", color=TEXT, fontsize=8)
    ax.set_ylabel("-log10 P (CellBender)", color=TEXT, fontsize=8)
    ax.set_title(title, color=TEXT, fontsize=10, loc="left", pad=22)
    ax.tick_params(colors=MUTED, labelsize=7)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    # legend above the plot so it never covers the most significant genes
    ax.legend(fontsize=7, frameon=False, loc="lower left", bbox_to_anchor=(0, 1.02), ncol=3, borderaxespad=0)
    fig.tight_layout()
    fig.savefig(path, facecolor=SURFACE)
    plt.close(fig)


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("cellbender_de_dir")
    p.add_argument("cellranger_de_dir")
    p.add_argument("--out", required=True)
    p.add_argument("--alpha", type=float, default=0.05)
    args = p.parse_args(argv)

    cb_dir, cr_dir, out = Path(args.cellbender_de_dir), Path(args.cellranger_de_dir), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for cb_file in sorted(cb_dir.glob("*_de.csv")):
        cr_file = cr_dir / cb_file.name
        if not cr_file.exists():
            print(f"skipping {cb_file.name}: no CellRanger result")
            continue
        label = cb_file.name.removesuffix("_de.csv")
        df = concordant(pd.read_csv(cb_file, index_col=0), pd.read_csv(cr_file, index_col=0), args.alpha)
        df.to_csv(out / f"{label}_concordant_de.csv")
        volcano(df, label, out / f"{label}_volcano.png")
        n_up, n_down = (df["direction"] == "up").sum(), (df["direction"] == "down").sum()
        rows.append({"celltype": label, "tested": len(df), "up": n_up, "down": n_down, "total_DE": n_up + n_down})
        print(f"{label:25s} DE in both = {n_up + n_down:5d}  (up {n_up}, down {n_down})")
    pd.DataFrame(rows).to_csv(out / "summary.csv", index=False)


if __name__ == "__main__":
    main()
