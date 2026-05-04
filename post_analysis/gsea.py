import pandas as pd
import numpy as np
from pathlib import Path
import gseapy as gp

from utils import reset_folder


def run_gsea(
    summary_csv,
    output_dir,
    celltype_key="celltype_pred",
    min_genes=10,
    permutation_num=1000,
    fdr_threshold=0.05,
    nes_threshold=1.5,
    top_n=5
):

    output_dir = Path(output_dir)
    reset_folder(output_dir)

    df = pd.read_csv(summary_csv)

    gene_sets = {
        "hallmark": "MSigDB_Hallmark_2020",
        "reactome": "Reactome_2022",
    }

    results_all = []

    # =========================
    # GSEA LOOP
    # =========================
    for ct in df[celltype_key].unique():

        print(f"\n=== CELLTYPE: {ct} ===")

        df_ct = df[df[celltype_key] == ct]
        ct_dir = output_dir / ct
        ct_dir.mkdir(exist_ok=True)

        for contrast in df_ct["contrast"].unique():

            print(f"  -> {contrast}")

            sub = df_ct[df_ct["contrast"] == contrast].copy()

            # -------------------------
            # CLEANING
            # -------------------------
            sub = sub.dropna(subset=["gene", "logfc", "fdr", "t_stat"])
            # REMOVE HOUSEKEEPING GENES
            sub = sub[~sub["gene"].str.match(r"^RPL|^RPS|^MT-")]

            # collapse duplicates
            sub = sub.groupby("gene", as_index=False).agg({
                "logfc": "mean",
                "fdr": "mean",
                "t_stat": "mean"
            })

            # =========================
            # RANKING (FIXED)
            # =========================
            sub["score"] = sub["t_stat"]

            # deterministic tie-breaking (no randomness)
            sub = sub.sort_values(
                ["score", "gene"],
                ascending=[False, True]
            )

            # enforce strictly unique ranking
            sub["score"] = sub["score"] + np.arange(len(sub)) * 1e-12

            if len(sub) < min_genes:
                continue

            rnk = sub[["gene", "score"]]

            for gs_name, gs_db in gene_sets.items():

                try:
                    enr = gp.prerank(
                        rnk=rnk,
                        gene_sets=gs_db,
                        processes=4,
                        min_size=10,
                        max_size=500,
                        permutation_num=permutation_num,
                        outdir=None,
                        seed=42,
                        verbose=False,
                    )

                    res = enr.res2d

                    if res is None or res.empty:
                        continue

                    res["celltype"] = ct
                    res["contrast"] = contrast
                    res["geneset"] = gs_name

                    results_all.append(res)

                    # save per contrast
                    res.to_csv(
                        ct_dir / f"{contrast}_{gs_name}.csv",
                        index=False
                    )

                except Exception as e:
                    print(f"ERROR ({ct}, {contrast}, {gs_name}):", e)
                    continue

    if len(results_all) == 0:
        raise RuntimeError("No enrichment results")

    # =========================
    # FINAL OUTPUT
    # =========================
    final = pd.concat(results_all, ignore_index=True)
    final.to_csv(output_dir / "ALL_ENRICHMENTS.csv", index=False)

    # =========================
    # FILTER SIGNIFICANT
    # =========================
    sig = final[
        (final["FDR q-val"] < fdr_threshold) &
        (final["NES"].abs() > nes_threshold)
    ].copy()

    # =========================
    # TOP PATHWAYS
    # =========================
    top = (
        sig.sort_values("NES", ascending=False)
           .groupby(["celltype", "contrast", "geneset"])
           .head(top_n)
    )

    top.to_csv(output_dir / "TOP_PATHWAYS.csv", index=False)

    return final, top