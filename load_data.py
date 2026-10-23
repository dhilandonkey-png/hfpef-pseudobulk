"""Load single-nucleus data into an AnnData with raw counts in .X.

Supports:
  * .h5ad files (e.g. an AnnData export of SCP3342)
  * Broad Single Cell Portal layout: a 10x-style matrix.mtx(.gz) with
    barcodes/features files, plus the SCP metadata.txt whose second line is
    the "TYPE" row (group / numeric).
"""

from __future__ import annotations

import gzip
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


def _open(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def read_scp_metadata(path: str | Path) -> pd.DataFrame:
    """Read an SCP metadata file, dropping the TYPE row if present."""
    path = Path(path)
    with _open(path) as fh:
        df = pd.read_csv(fh, sep="\t", dtype=str)
    if len(df) and str(df.iloc[0, 0]).upper() == "TYPE":
        types = df.iloc[0]
        df = df.iloc[1:].copy()
        for col, kind in types.items():
            if str(kind).lower() == "numeric":
                df[col] = pd.to_numeric(df[col], errors="coerce")
    return df.set_index(df.columns[0])


def read_mtx_dir(matrix: Path, barcodes: Path, features: Path) -> ad.AnnData:
    with _open(matrix) as fh:
        mat = scipy.io.mmread(fh)
    mat = sp.csr_matrix(mat)
    with _open(barcodes) as fh:
        bc = [line.strip().split("\t")[0] for line in fh if line.strip()]
    with _open(features) as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    # 10x features.tsv: gene_id, gene_symbol, feature_type -> prefer symbol
    genes = [r[1] if len(r) > 1 else r[0] for r in rows]

    # 10x / SCP matrices are genes x cells; transpose to cells x genes
    if mat.shape == (len(genes), len(bc)):
        mat = mat.T.tocsr()
    elif mat.shape != (len(bc), len(genes)):
        raise ValueError(
            f"matrix shape {mat.shape} matches neither "
            f"{len(genes)} genes x {len(bc)} barcodes nor the transpose"
        )
    adata = ad.AnnData(X=mat, obs=pd.DataFrame(index=bc), var=pd.DataFrame(index=genes))
    adata.var_names_make_unique()
    return adata


def read_h5ad_layer(path: str, layer: str) -> ad.AnnData:
    """Read obs, var and a single layer as X, skipping X and the other layers.

    SCP3342's h5ad holds normalised X plus two raw-count layers; loading only
    the one we need keeps memory to roughly a third.
    """
    import h5py
    from anndata.io import read_elem

    with h5py.File(path, "r") as f:
        if layer not in f["layers"]:
            raise ValueError(f"layer '{layer}' not found; available: {list(f['layers'].keys())}")
        return ad.AnnData(X=read_elem(f["layers"][layer]), obs=read_elem(f["obs"]), var=read_elem(f["var"]))


def load(
    h5ad: str | None = None,
    matrix: str | None = None,
    barcodes: str | None = None,
    features: str | None = None,
    metadata: str | None = None,
    layer: str | None = None,
) -> ad.AnnData:
    if h5ad and layer:
        adata = read_h5ad_layer(h5ad, layer)
    elif h5ad:
        adata = ad.read_h5ad(h5ad)
        if adata.raw is not None and not _looks_like_counts(adata.X):
            # Many published objects keep log-normalised data in X and counts in .raw
            adata = ad.AnnData(X=adata.raw.X, obs=adata.obs, var=adata.raw.var)
    elif matrix and barcodes and features:
        adata = read_mtx_dir(Path(matrix), Path(barcodes), Path(features))
    else:
        raise ValueError("provide --h5ad, or --matrix/--barcodes/--features")

    if metadata:
        meta = read_scp_metadata(metadata)
        shared = adata.obs_names.intersection(meta.index)
        if len(shared) == 0:
            raise ValueError("no barcodes shared between matrix and metadata")
        adata = adata[shared].copy()
        adata.obs = adata.obs.join(meta, how="left", rsuffix="_meta")

    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    return adata


def _looks_like_counts(X, n: int = 10000) -> bool:
    vals = X.data if sp.issparse(X) else np.asarray(X).ravel()
    vals = vals[:n]
    return bool(len(vals) == 0 or (np.all(vals >= 0) and np.allclose(vals, np.round(vals))))
