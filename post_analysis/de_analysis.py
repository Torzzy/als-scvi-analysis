import os
import tempfile
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .io import load_adata
from .utils import reset_folder


# ======================================================
# R SCRIPT (ROBUST EDGE R DESIGN + GSEA RANKING)
# ======================================================
R_SCRIPT = r"""
suppressMessages(library(edgeR))

meta <- read.csv("__META__")
counts <- read.csv("__COUNTS__", row.names=1)

meta$condition <- factor(meta$condition)
meta$region <- factor(meta$region)

design <- model.matrix(~ condition * region, data=meta)

y <- DGEList(counts=counts)

keep <- filterByExpr(y, design)
y <- y[keep,,keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design)

# ======================================================
# SAFE TEST FUNCTION (WITH GSEA RANKING)
# ======================================================
safe_test <- function(coef_name, out_name) {

    if (!(coef_name %in% colnames(design))) {
        cat("SKIP:", coef_name, "\n")
        return(NULL)
    }

    res <- glmQLFTest(fit, coef=coef_name)
    tt <- topTags(res, n=Inf)$table
    tt$gene <- rownames(tt)

    # FORCE FDR ALWAYS
    tt$FDR <- p.adjust(tt$PValue, method="BH")

    # ======================================================
    # GSEA RANKING FEATURES (ADDED ONLY, NOTHING BROKEN)
    # ======================================================

    tt$logFC_rank <- tt$logFC

    # pseudo t-statistic (edgeR standard proxy)
    if ("F" %in% colnames(tt)) {
        tt$t_stat <- sign(tt$logFC) * sqrt(tt$F)
    } else if ("LR" %in% colnames(tt)) {
        tt$t_stat <- sign(tt$logFC) * sqrt(tt$LR)
    } else {
        tt$t_stat <- tt$logFC
    }

    # FINAL GSEA RANKING SCORE
    tt$rank_score <- tt$t_stat

    write.csv(tt, paste0("__OUT__", "_", out_name, ".csv"), row.names=FALSE)
}

# ======================================================
# CONTRASTS
# ======================================================

# 1) GLOBAL ALS vs CTRL
safe_test("conditionCTRL", "global_condition")

# 2) REGION EFFECTS
for (r in c("FX","MCX","SC")) {
    safe_test(paste0("region", r), paste0("region_", r))
}

# 3) INTERACTIONS ALS × REGION
for (r in c("FX","MCX","SC")) {
    safe_test(
        paste0("conditionCTRL:region", r),
        paste0("interaction_", r)
    )
}
"""

# ======================================================
# MAIN PIPELINE
# ======================================================
def run_pseudobulk_de(
    latent_path,
    output_dir,
    celltype_key="celltype_pred",
    condition_key="condition",
    region_key="region",
    patient_key="patient_id",
    group1="ALS",
    group2="CTRL",
    min_cells=15,
    min_patients=2,
):

    output_dir = Path(output_dir)
    reset_folder(output_dir)

    adata = load_adata(latent_path)

    obs = adata.obs
    genes = adata.var_names.to_numpy()

    X = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)

    for ct in obs[celltype_key].unique():

        print(f"\n=== {ct} ===")

        ct_dir = output_dir / str(ct)
        ct_dir.mkdir(parents=True, exist_ok=True)

        mask = (obs[celltype_key] == ct).to_numpy()
        idx_ct = np.where(mask)[0]

        if len(idx_ct) == 0:
            continue

        patient = obs[patient_key].to_numpy()[idx_ct]
        condition = obs[condition_key].to_numpy()[idx_ct]
        region = obs[region_key].to_numpy()[idx_ct]

        patient_vectors = {}
        meta_rows = []

        for p in np.unique(patient):

            p_mask = (patient == p)
            cell_idx = idx_ct[p_mask]

            if len(cell_idx) < min_cells:
                continue

            vec = np.asarray(X[cell_idx].sum(axis=0)).ravel()

            patient_vectors[p] = vec

            meta_rows.append({
                "sample": p,
                "condition": condition[p_mask][0],
                "region": region[p_mask][0],
            })

        if len(patient_vectors) < 4:
            continue

        meta = pd.DataFrame(meta_rows)

        if (meta["condition"] == group1).sum() < min_patients:
            continue
        if (meta["condition"] == group2).sum() < min_patients:
            continue

        counts = np.vstack(list(patient_vectors.values()))

        with tempfile.TemporaryDirectory() as tmp:

            counts_f = os.path.join(tmp, "counts.csv")
            meta_f = os.path.join(tmp, "meta.csv")
            out_prefix = os.path.join(tmp, "out")
            script_f = os.path.join(tmp, "run.R")

            pd.DataFrame(
                counts.T,
                index=genes,
                columns=meta["sample"]
            ).to_csv(counts_f)

            meta.to_csv(meta_f, index=False)

            script = R_SCRIPT.replace("__META__", meta_f)
            script = script.replace("__COUNTS__", counts_f)
            script = script.replace("__OUT__", out_prefix)

            with open(script_f, "w") as f:
                f.write(script)

            res = subprocess.run(
                ["Rscript", script_f],
                capture_output=True,
                text=True
            )

            print(res.stdout)
            print(res.stderr)

            if res.returncode != 0:
                raise RuntimeError(f"R failed for {ct}")

            for f in os.listdir(tmp):
                if f.endswith(".csv"):
                    df = pd.read_csv(os.path.join(tmp, f))
                    df.to_csv(ct_dir / f, index=False)
# ======================================================
# SUMMARY FUNCTION
# ======================================================
from pathlib import Path
import pandas as pd


def summarize_edger_results(
    de_dir,
    output_csv,
    celltype_key="celltype_pred"
):

    de_dir = Path(de_dir)

    files = list(de_dir.rglob("out_*.csv"))

    print("FILES FOUND:", files)

    all_dfs = []

    for f in files:

        if "meta" in f.name or "counts" in f.name:
            continue

        df = pd.read_csv(f)

        if "gene" not in df.columns:
            continue

        celltype = f.parent.name
        contrast = f.stem.replace("out_", "")

        df.columns = [c.lower() for c in df.columns]

        if "fdr" not in df.columns:
            if "padj" in df.columns:
                df["fdr"] = df["padj"]
            else:
                df["fdr"] = 1.0

        # ==============================
        # ADD t-statistic (IF PRESENT)
        # ==============================
        if "t_stat" not in df.columns:
            if "rank_score" in df.columns:
                df["t_stat"] = df["rank_score"]
            else:
                df["t_stat"] = df["logfc"]

        sub = df[["gene", "logfc", "fdr", "t_stat"]].copy()
        sub[celltype_key] = celltype
        sub["contrast"] = contrast

        all_dfs.append(sub)

    if len(all_dfs) == 0:
        raise RuntimeError("No valid results found")

    merged = pd.concat(all_dfs, axis=0, ignore_index=True)

    merged.to_csv(output_csv, index=False)

    return merged