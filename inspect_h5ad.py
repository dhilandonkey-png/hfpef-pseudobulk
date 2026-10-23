"""Print what is inside an h5ad so you can pick the pseudobulk column names.

  python -I inspect_h5ad.py SCP3342.h5ad
"""

from __future__ import annotations

import sys

import anndata as ad
import numpy as np
import scipy.sparse as sp


def looks_like_counts(X, n: int = 100000) -> bool:
    vals = X.data[:n] if sp.issparse(X) else np.asarray(X[:200]).ravel()[:n]
    return bool(len(vals) == 0 or (np.all(vals >= 0) and np.allclose(vals, np.round(vals))))


def main(path: str) -> None:
    adata = ad.read_h5ad(path, backed="r")
    print(f"{adata.n_obs:,} nuclei x {adata.n_vars:,} genes\n")

    print("obs columns (candidates for --donor-col / --celltype-col / --group-col / --covariates):")
    for col in adata.obs.columns:
        s = adata.obs[col]
        n = s.nunique()
        preview = ", ".join(map(str, s.astype(str).value_counts().index[:6]))
        print(f"  {col:30s} {str(s.dtype):10s} {n:7,d} unique   e.g. {preview}")

    X = adata.X[:2000]
    print(f"\nX: {type(X).__name__}, raw integer counts: {looks_like_counts(X)}")
    for name in (k for k in adata.layers.keys() if k is not None):  # backed mode can list None
        print(f"layer '{name}': raw integer counts: {looks_like_counts(adata.layers[name][:2000])}")
    if adata.raw is not None:
        print(f".raw: {adata.raw.shape}, raw integer counts: {looks_like_counts(adata.raw.X[:2000])}")
    print(f"obsm: {list(adata.obsm.keys())}")


if __name__ == "__main__":
    main(sys.argv[1])
