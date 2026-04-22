import os
import glob
import numpy as np
import pandas as pd
import gseapy as gp

import utils


def _ranking(df):
    """
    Old behavior:
    score = sign(logFC) * -log10(fdr)
    """

    d = df.copy()

    if "gene" not in d.columns:
        if "index" in d.columns:
            d = d.rename(columns={"index": "gene"})
        else:
            d = d.reset_index().rename(columns={"index": "gene"})

    required = ["gene", "logFC", "fdr"]
    for col in required:
        if col not in d.columns:
            raise ValueError(f"Missing required column: {col}")

    d = d.dropna(subset=required).copy()

    d["gene"] = d["gene"].astype(str)
    d["logFC"] = pd.to_numeric(d["logFC"], errors="coerce")
    d["fdr"] = pd.to_numeric(d["fdr"], errors="coerce")
    d = d.dropna(subset=["logFC", "fdr"])

    d["fdr"] = d["fdr"].clip(lower=1e-300, upper=1.0)

    d["score"] = np.sign(d["logFC"]) * (-np.log10(d["fdr"]))
    d = d[np.isfinite(d["score"])]

    d = d.drop_duplicates(subset="gene", keep="first")

    rng = np.random.default_rng(42)
    d["score"] += rng.normal(0, 1e-9, len(d))

    d = d.sort_values("score", ascending=False)

    return d[["gene", "score"]]


def _sort_results(df):
    if "FDR q-val" in df.columns:
        return df.sort_values(["FDR q-val", "NES"], ascending=[True, False])
    if "NES" in df.columns:
        return df.sort_values("NES", ascending=False)
    return df


def _top_pathways(df, n=10):
    if df.empty:
        return pd.DataFrame()

    d = _sort_results(df.copy())

    keep = [
        c for c in [
            "comparison",
            "cell_type",
            "region",
            "geneset_db",
            "Term",
            "NES",
            "FDR q-val",
        ]
        if c in d.columns
    ]

    d = d[keep].head(n).copy()

    d = d.rename(columns={
        "geneset_db": "database",
        "Term": "pathway",
        "FDR q-val": "FDR",
    })

    if "NES" in d.columns:
        d["direction"] = np.where(d["NES"] > 0, "ALS_up", "CTRL_up")

    return d


def _parse_comparison_name(name):
    """
    Parse:
    Excitatory_FCX_DE -> (Excitatory, FCX)
    Astrocyte_Motor_Cortex_DE -> (Astrocyte, Motor_Cortex)
    """

    stem = name.replace("_DE", "")

    known_regions = [
        "FCX",
        "TCX",
        "SC",
        "HIP",
        "FX",
        "MCX",
        "Motor_Cortex",
        "Spinal_Cord",
        "Frontal_Cortex",
        "Temporal_Cortex",
    ]

    for rg in sorted(known_regions, key=len, reverse=True):
        suffix = "_" + rg
        if stem.endswith(suffix):
            return stem[:-len(suffix)], rg

    parts = stem.split("_")
    if len(parts) >= 2:
        return "_".join(parts[:-1]), parts[-1]

    return stem, "Unknown"


def run_gsea(
    de_dir,
    output_dir,
    gene_sets=("MSigDB_Hallmark_2020", "Reactome_2022"),
    min_genes_ranked=100,
    min_sig_genes=0,
    fdr_threshold=0.25,
    min_size=10,
    max_size=500,
    permutation_num=1000,
    processes=4,
    top_n=10,
):
    """
    Exécute une analyse GSEA (Gene Set Enrichment Analysis) à partir de
    résultats d'expression différentielle.

    La fonction parcourt les fichiers DE, construit un ranking de gènes,
    lance un preranked GSEA pour chaque base de gènes fournie, filtre les
    voies significatives, puis sauvegarde les résultats détaillés ainsi que
    des tableaux résumés des top pathways enrichis.

    :param de_dir: Dossier contenant les fichiers de résultats DE
        (`*_DE.csv`).
    :type de_dir: str | Path

    :param output_dir: Dossier de sortie des résultats GSEA.
    :type output_dir: str | Path

    :param gene_sets: Bases de signatures biologiques utilisées pour GSEA.
    :type gene_sets: tuple[str, ...]

    :param min_genes_ranked: Nombre minimal de gènes requis dans le ranking
        pour lancer GSEA.
    :type min_genes_ranked: int

    :param min_sig_genes: Nombre minimal de gènes significatifs (selon FDR)
        requis pour analyser une comparaison.
    :type min_sig_genes: int

    :param fdr_threshold: Seuil de significativité utilisé pour filtrer
        les résultats GSEA.
    :type fdr_threshold: float

    :param min_size: Taille minimale des gene sets testés.
    :type min_size: int

    :param max_size: Taille maximale des gene sets testés.
    :type max_size: int

    :param permutation_num: Nombre de permutations utilisées par GSEA.
    :type permutation_num: int

    :param processes: Nombre de processus parallèles utilisés.
    :type processes: int

    :param top_n: Nombre de pathways principaux à conserver par comparaison.
    :type top_n: int

    :return: Tableau récapitulatif de l'ensemble des enrichissements détectés.
    :rtype: pandas.DataFrame
    """

    utils.reset_folder(output_dir)

    files = glob.glob(os.path.join(de_dir, "*_DE.csv"))
    if not files:
        print("No DE files found.")
        return pd.DataFrame()

    all_results = []

    for fp in files:

        name = os.path.basename(fp).replace(".csv", "")
        print(f"\n=== {name} ===")

        try:
            df = pd.read_csv(fp)
        except Exception as e:
            print(f"Read failed: {e}")
            continue

        if len(df) < min_genes_ranked:
            print("Skipped: too few genes")
            continue

        if "fdr" not in df.columns:
            print("Skipped: no fdr column")
            continue

        n_sig = int((pd.to_numeric(df["fdr"], errors="coerce") < fdr_threshold).sum())

        if n_sig < min_sig_genes:
            print(f"Skipped: only {n_sig} significant genes")
            continue

        try:
            rnk = _ranking(df)
        except Exception as e:
            print(f"Ranking failed: {e}")
            continue

        if len(rnk) < min_genes_ranked:
            print("Skipped: ranking too small")
            continue

        cell_type, region = _parse_comparison_name(name)
        comp_results = []

        for gs in gene_sets:

            print(f"Running {gs}")

            try:
                pre = gp.prerank(
                    rnk=rnk,
                    gene_sets=gs,
                    permutation_num=permutation_num,
                    min_size=min_size,
                    max_size=max_size,
                    processes=processes,
                    seed=42,
                    outdir=None,
                    verbose=False,
                )

                res = pre.res2d

                if res is None or res.empty:
                    continue

                res = res.copy()

                if "Term" not in res.columns:
                    res["Term"] = res.index.astype(str)

                if "FDR q-val" not in res.columns or "NES" not in res.columns:
                    continue

                # old robust filtering
                res = res[
                    (res["FDR q-val"] < fdr_threshold) &
                    (res["NES"].abs() > 1.2)
                ].copy()

                if res.empty:
                    print("No robust pathway")
                    continue

                res["comparison"] = name
                res["cell_type"] = cell_type
                res["region"] = region
                res["geneset_db"] = gs

                res = _sort_results(res)

                comp_results.append(res)
                all_results.append(res)

                out_fp = os.path.join(
                    output_dir,
                    f"{name}_{gs}_GSEA.csv"
                )
                res.to_csv(out_fp, index=False)

            except Exception as e:
                print(f"Failed {gs}: {e}")

        if comp_results:
            comp_df = pd.concat(comp_results, ignore_index=True)
            top = _top_pathways(comp_df, n=top_n)

            top_fp = os.path.join(
                output_dir,
                f"{name}_TOP_pathways.csv"
            )
            top.to_csv(top_fp, index=False)

    if not all_results:
        print("No GSEA results generated.")
        return pd.DataFrame()

    summary = pd.concat(all_results, ignore_index=True)
    summary = _sort_results(summary)

    full_fp = os.path.join(output_dir, "GSEA_summary_results.csv")
    summary.to_csv(full_fp, index=False)

    # one readable global table
    top_rows = []

    for (comp, db), sub in summary.groupby(["comparison", "geneset_db"]):

        sub = _sort_results(sub).head(top_n).copy()

        for _, row in sub.iterrows():
            top_rows.append({
                "comparison": comp,
                "cell_type": row["cell_type"],
                "region": row["region"],
                "database": db,
                "pathway": row["Term"],
                "NES": row["NES"],
                "FDR": row["FDR q-val"],
                "direction": "ALS_up" if row["NES"] > 0 else "CTRL_up",
            })

    top_global = pd.DataFrame(top_rows)

    if not top_global.empty:
        top_global = top_global.sort_values(
            ["FDR", "NES"],
            ascending=[True, False]
        )

        top_fp = os.path.join(output_dir, "GSEA_TOP_pathways.csv")
        top_global.to_csv(top_fp, index=False)

        print("\n========== TOP PATHWAYS ==========\n")
        print(top_global.to_string(index=False))

    print("\nSaved:")
    print(full_fp)

    return summary