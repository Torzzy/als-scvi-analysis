import scipy.sparse as sp
import subprocess
from pathlib import Path
import os
import numpy as np
import pandas as pd
from tempfile import TemporaryDirectory


R_SCRIPT = r'''
suppressMessages(library(edgeR))

args <- commandArgs(trailingOnly=TRUE)

counts <- read.csv(args[1], row.names=1, check.names=FALSE)
meta   <- read.csv(args[2], stringsAsFactors=FALSE)
outf   <- args[3]
g1     <- args[4]
g2     <- args[5]

# ----------------------------
# sanity checks
# ----------------------------
stopifnot(all(c("sample", "group") %in% colnames(meta)))

meta$group <- factor(meta$group, levels=c(g2, g1))

# ----------------------------
# edgeR pipeline
# ----------------------------
y <- DGEList(counts=as.matrix(counts))

# filter lowly expressed genes
keep <- filterByExpr(y, group=meta$group)
y <- y[keep, , keep.lib.sizes=FALSE]

# normalization
y <- calcNormFactors(y)

# design matrix
design <- model.matrix(~ group, data=meta)

# dispersion + GLM
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

qlf <- glmQLFTest(fit, coef=2)

res <- topTags(qlf, n=Inf)$table

# add gene names
res$gene <- rownames(res)

# clean ordering
res <- res[order(res$FDR, -abs(res$logFC)), ]

write.csv(res, outf, row.names=FALSE)
'''

def build_pseudobulk_matrix(adata, groupby=("cell_type", "region"), patient_key="patient_id"):
    """
    Construit une matrice pseudobulk en agrégeant les comptes d’expression
    au niveau de groupes définis (ex: celltype, région) et par patient.

    Étapes :
    - Récupère la matrice d’expression (dense ou sparse)
    - Définit des groupes basés sur les colonnes spécifiées dans `groupby`
      ainsi que l’identifiant patient
    - Pour chaque combinaison (groupe, patient) :
        - somme les expressions des cellules correspondantes
    - Construit une matrice pseudobulk (échantillons x gènes)
      et une table de métadonnées associée

    :param adata: objet AnnData contenant les données (cells x genes)
    :param groupby: tuple de colonnes dans adata.obs pour définir les groupes
    :param patient_key: colonne dans adata.obs identifiant les patients
    :return:
        - X_pb: matrice pseudobulk (n_groups x n_genes)
        - meta: DataFrame des métadonnées (groupes + patient)
        - genes: noms des gènes
    """

    obs = adata.obs
    X = adata.X.tocsr() if sp.issparse(adata.X) else adata.X
    genes = adata.var_names.to_numpy()

    keys = list(zip(*(obs[c].to_numpy() for c in groupby)))
    patients = obs[patient_key].to_numpy()

    df = pd.DataFrame({
        "key": keys,
        "patient": patients
    })

    groups = df.groupby(["key", "patient"]).indices

    rows = []
    meta = []

    for (key, patient), idx in groups.items():

        vec = np.asarray(X[idx].sum(axis=0)).ravel()

        rows.append(vec)
        meta.append((*key, patient))

    meta = pd.DataFrame(meta, columns=list(groupby) + ["patient"])
    X_pb = np.vstack(rows)

    return X_pb, meta, genes

def run_cluster_vs_rest_de_fast(X_pb, meta, genes, output_dir, cluster_key="cell_type"):
    """
    Lance une analyse différentielle (cluster vs reste) à partir d’une matrice pseudobulk.

    :param X_pb: matrice pseudobulk (samples x genes ou gènes x samples selon usage)
    :param meta: DataFrame contenant les métadonnées des pseudobulks
    :param genes: liste des noms de gènes
    :param output_dir: dossier de sortie pour les résultats
    :param cluster_key: colonne de meta définissant les clusters
    :return: aucun retour (sauvegarde des résultats sur disque)
    """

    clusters = meta[cluster_key].unique()

    os.makedirs(output_dir, exist_ok=True)

    for cl in clusters:

        mask = meta[cluster_key].to_numpy() == cl

        group = np.where(mask, "cluster", "rest")

        if mask.sum() < 5:
            continue

        with TemporaryDirectory() as tmp:

            cf = os.path.join(tmp, "counts.csv")
            mf = os.path.join(tmp, "meta.csv")
            rf = os.path.join(tmp, "run.R")
            of = os.path.join(tmp, "out.csv")

            pd.DataFrame(X_pb.T, index=genes).to_csv(cf)

            pd.DataFrame({
                "sample": np.arange(len(meta)),
                "group": group
            }).to_csv(mf, index=False)

            with open(rf, "w") as f:
                f.write(R_SCRIPT)

            subprocess.run(
                ["Rscript", rf, cf, mf, of, "cluster", "rest"],
                check=True
            )

            if os.path.exists(of):
                res = pd.read_csv(of)
                res.to_csv(os.path.join(output_dir, f"{cl}_vs_rest.csv"), index=False)


def annotate_clusters_from_de(
    de_folder,
    marker_dict,
    fdr_col="FDR",
    logfc_col="logFC",
    eps=1e-10,
    margin=0.2,
):
    """
    Annoter des clusters à partir de résultats d’expression différentielle
    (cluster vs reste), en utilisant des gènes marqueurs.

    Étapes :
    - Parcourt les fichiers DE dans un dossier (format *_vs_rest.csv)
    - Pour chaque cluster :
        - vérifie la présence des colonnes nécessaires
        - construit un score de gène basé sur logFC et significativité (FDR)
        - agrège ces scores par type cellulaire via les gènes marqueurs
    - Détermine le type cellulaire dominant :
        - basé sur le meilleur score moyen des marqueurs
        - avec une règle de séparation (margin) pour éviter les assignations ambiguës
    - Calcule une mesure de confiance basée sur le ratio du meilleur score

    :param de_folder: dossier contenant les résultats DE par cluster
    :param marker_dict: dict {celltype -> liste de gènes marqueurs}
    :param fdr_col: nom de la colonne FDR dans les fichiers DE
    :param logfc_col: nom de la colonne log fold-change
    :param eps: petite constante pour éviter les divisions/log de zéro
    :param margin: seuil de différence entre meilleur et second score pour éviter les ambiguïtés
    :return: DataFrame avec annotation des clusters
    """

    de_folder = Path(de_folder)

    files = sorted(de_folder.glob("*_vs_rest.csv"))

    if len(files) == 0:
        raise ValueError(f"No DE files found in {de_folder}")

    results = []

    for f in files:
        df = pd.read_csv(f)
        required = {"gene", logfc_col, fdr_col}
        if not required.issubset(df.columns):
            print(f"[SKIP] {f.name} missing {required - set(df.columns)}")
            continue

        df = df.copy()

        df[fdr_col] = df[fdr_col].replace(0, eps)

        df["score"] = df[logfc_col] * (-np.log10(df[fdr_col] + eps))
        gene_score = dict(zip(df["gene"], df["score"]))

        type_scores = {}
        marker_coverage = {}

        for cell_type, markers in marker_dict.items():

            scores = [gene_score[g] for g in markers if g in gene_score]

            type_scores[cell_type] = float(np.mean(scores)) if scores else 0.0
            marker_coverage[cell_type] = len(scores)

        if not type_scores:
            continue

        sorted_types = sorted(type_scores.items(), key=lambda x: x[1], reverse=True)

        best_type, best_score = sorted_types[0]
        second_score = sorted_types[1][1] if len(sorted_types) > 1 else 0.0

        score_diff = best_score - second_score

        if best_score < 1e-6:
            label = "NoSignal"
        elif score_diff < margin:
            label = "Mixed"
        else:
            label = best_type

        # confidence metric (IMPORTANT)
        confidence = best_score / (second_score + eps + best_score)

        cluster_name = f.stem.replace("_vs_rest", "")

        results.append({
            "cluster": cluster_name,
            "assigned_type": label,
            "best_type": best_type,
            "best_score": best_score,
            "second_score": second_score,
            "margin": score_diff,
            "confidence": confidence,
            "top_marker_coverage": marker_coverage[best_type] / max(len(marker_dict[best_type]), 1)
        })

    res_df = pd.DataFrame(results)

    if res_df.empty:
        raise RuntimeError(
            f"No valid DE results parsed in {de_folder}. "
            "Check that CSVs contain gene/logFC/FDR columns."
        )

    res_df = res_df.sort_values("confidence", ascending=False)

    return res_df