"""Generate a small fake dataset shaped like SCP3342 to test the pipeline.

19 HFpEF + 24 control donors in pools of 6, a handful of cell types with
realistic imbalance, and a known set of genes shifted in HFpEF cardiomyocytes.
This is NOT real data.
"""

from __future__ import annotations

import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

CELLTYPES = {  # cell type -> mean nuclei per donor
    "Cardiomyocyte": 300,
    "Fibroblast": 270,
    "Endothelial1": 190,
    "Pericyte": 120,
    "Macrophage": 60,
    "Lymphocyte": 15,  # too rare in most donors -> should fail the DE filter
}
N_GENES = 2000
TRUE_DE = 100  # first 100 genes are up 2x in HFpEF cardiomyocytes


def main(out_path: str, seed: int = 0) -> None:
    rng = np.random.default_rng(seed)
    donors = [f"HFpEF_{i:02d}" for i in range(19)] + [f"CTRL_{i:02d}" for i in range(24)]
    order = rng.permutation(len(donors))
    pool = {donors[j]: f"pool{k // 6 + 1}" for k, j in enumerate(order)}
    sex = {d: rng.choice(["F", "M"]) for d in donors}

    base = {ct: rng.gamma(0.3, 2.0, N_GENES) for ct in CELLTYPES}
    blocks, obs = [], []
    for d in donors:
        disease = "HFpEF" if d.startswith("HFpEF") else "control"
        donor_effect = rng.lognormal(0, 0.2, N_GENES)
        for ct, mean_n in CELLTYPES.items():
            n = rng.poisson(mean_n)
            if n == 0:
                continue
            mu = base[ct] * donor_effect
            if disease == "HFpEF" and ct == "Cardiomyocyte":
                mu = mu.copy()
                mu[:TRUE_DE] *= 2.0
            depth = rng.lognormal(0, 0.3, n)[:, None]
            blocks.append(sp.csr_matrix(rng.poisson(mu[None, :] * depth)))
            obs += [(d, ct, disease, sex[d], pool[d])] * n

    X = sp.vstack(blocks).tocsr()
    obs = pd.DataFrame(obs, columns=["donor_id", "cell_type", "disease", "sex", "pool"])
    obs.index = [f"cell{i}" for i in range(len(obs))]
    var = pd.DataFrame(index=[f"GENE{i:04d}" for i in range(N_GENES)])
    ad.AnnData(X=X.astype(np.float32), obs=obs, var=var).write_h5ad(out_path)
    print(f"wrote {X.shape[0]:,} synthetic nuclei to {out_path}")


if __name__ == "__main__":
    main(sys.argv[1])
