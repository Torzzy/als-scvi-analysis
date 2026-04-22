import os
import tempfile
import subprocess
import joblib
import numpy as np
import pandas as pd
import scipy.sparse as sp

from .io import load_adata
from .utils import reset_folder


R_SCRIPT = r'''
suppressMessages(library(edgeR))

args <- commandArgs(trailingOnly=TRUE)
counts <- read.csv(args[1], row.names=1, check.names=FALSE)
meta   <- read.csv(args[2], stringsAsFactors=FALSE)
outf   <- args[3]
g1     <- args[4]
g2     <- args[5]

meta$group <- factor(meta$group, levels=c(g2, g1))

y <- DGEList(counts=counts)
keep <- filterByExpr(y, group=meta$group)
y <- y[keep,,keep.lib.sizes=FALSE]

y <- calcNormFactors(y)

design <- model.matrix(~ group, data=meta)

y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
res <- topTags(glmQLFTest(fit, coef=2), n=Inf)$table

res$gene <- rownames(res)
write.csv(res, outf, row.names=FALSE)
'''


def run_pseudobulk_de(
    latent_path,
    classifier_path,
    output_dir,
    groupby="condition",
    split_by=("cell_type", "region"),
    group1="ALS",
    group2="CTRL",
    patient_key="patient_id",
    min_cells=15,
    min_patients=2,
    min_total_cells=30,
    conf_score=0.95,
    conf_margin=0.50,
    reclassify=True
):
    """
    Réalise une analyse d'expression différentielle pseudobulk entre deux groupes
    (ex. ALS vs CTRL) à partir de données single-cell annotées.

    Les cellules sont d'abord annotées via un classifieur entraîné sur l'espace
    latent scVI. Les cellules sont ensuite filtrées selon leur niveau de confiance,
    agrégées par patient pour former des profils pseudobulk, puis comparées avec
    edgeR pour chaque combinaison définie par `split_by`.

    :param latent_path: Chemin vers le fichier .h5ad contenant les données
        et la représentation latente scVI.
    :type latent_path: str | Path

    :param classifier_path: Chemin vers le classifieur sauvegardé (.pkl).
    :type classifier_path: str | Path

    :param output_dir: Dossier de sortie des résultats DE.
    :type output_dir: str | Path

    :param groupby: Colonne de `adata.obs` définissant les groupes à comparer.
    :type groupby: str

    :param split_by: Colonnes utilisées pour stratifier les analyses
        (ex. type cellulaire, région).
    :type split_by: tuple[str, ...]

    :param group1: Nom du premier groupe de comparaison.
    :type group1: str

    :param group2: Nom du second groupe de comparaison.
    :type group2: str

    :param patient_key: Colonne identifiant les patients / échantillons.
    :type patient_key: str

    :param min_cells: Nombre minimal de cellules par patient pour construire
        un pseudobulk.
    :type min_cells: int

    :param min_patients: Nombre minimal de patients par groupe requis
        pour effectuer le test statistique.
    :type min_patients: int

    :param min_total_cells: Nombre minimal total de cellules requis
        pour analyser un sous-ensemble.
    :type min_total_cells: int

    :param conf_score: Score minimal de probabilité du classifieur
        pour conserver une cellule.
    :type conf_score: float

    :param conf_margin: Écart minimal entre la meilleure et la seconde
        probabilité prédite pour conserver une cellule.
    :type conf_margin: float

    :return: None
    :rtype: None
    """

    reset_folder(output_dir)

    adata = load_adata(latent_path)
    adata.obs_names_make_unique()

    clfobj = joblib.load(classifier_path)
    clf = clfobj["model"]
    le = clfobj.get("label_encoder", None)

    X_latent = adata.obsm["X_scVI"]

    pred_enc = clf.predict(X_latent)
    pred = le.inverse_transform(pred_enc) if le is not None else pred_enc
    adata.obs["cell_type"] = pred

    # --------------------------------------------------
    # Reclassification + confidence mask
    # --------------------------------------------------
    if reclassify:
        clfobj = joblib.load(classifier_path)
        clf = clfobj["model"]
        le = clfobj.get("label_encoder", None)

        X_latent = adata.obsm["X_scVI"]

        pred_enc = clf.predict(X_latent)
        pred = le.inverse_transform(pred_enc) if le is not None else pred_enc
        adata.obs["cell_type"] = pred

        if hasattr(clf, "predict_proba"):
            probs = clf.predict_proba(X_latent)

            top1_idx = np.argmax(probs, axis=1)
            top1_score = probs[np.arange(len(probs)), top1_idx]

            probs_sorted = np.sort(probs, axis=1)
            top2_score = probs_sorted[:, -2]

            margin = top1_score - top2_score
            conf_mask = (
                    (top1_score > conf_score) &
                    (margin > conf_margin)
            )
        else:
            conf_mask = np.ones(adata.n_obs, dtype=bool)

    else:
        # keep all cells, no reclassification
        conf_mask = np.ones(adata.n_obs, dtype=bool)

    genes = adata.var_names.to_numpy()

    # --------------------------------------------------
    # Main loop
    # --------------------------------------------------
    grouped = adata.obs.groupby(list(split_by), observed=False).indices

    for keys, obs_idx in grouped.items():

        idx = np.asarray(obs_idx, dtype=int)

        idx = idx[conf_mask[idx]]

        # old behavior: skip small subsets
        if len(idx) < min_total_cells:
            continue

        sub = adata[idx].copy()
        obs = sub.obs.copy()

        # ----------------------------------------------
        # keep patients with enough cells
        # ----------------------------------------------
        pats = obs[patient_key].value_counts()
        keep_pat = pats[pats >= min_cells].index

        keep_mask = obs[patient_key].isin(keep_pat).values
        sub = sub[keep_mask].copy()
        obs = sub.obs.copy()

        if sub.n_obs == 0:
            continue

        # ----------------------------------------------
        # build pseudobulk
        # ----------------------------------------------
        rows = []
        meta_rows = []

        for p in obs[patient_key].unique():

            pm = (obs[patient_key] == p).values
            Xp = sub.X[pm]

            if Xp.shape[0] < min_cells:
                continue

            if sp.issparse(Xp):
                vec = np.asarray(Xp.sum(axis=0)).ravel()
            else:
                vec = np.asarray(Xp.sum(axis=0)).ravel()

            rows.append(vec)

            row0 = obs.loc[pm].iloc[0]
            meta_rows.append(
                {
                    "sample": p,
                    "group": row0[groupby],
                }
            )

        if len(rows) == 0:
            continue

        meta = pd.DataFrame(meta_rows)

        g1 = (meta["group"] == group1).sum()
        g2 = (meta["group"] == group2).sum()

        if g1 < min_patients or g2 < min_patients:
            continue

        counts = np.vstack(rows)

        # ----------------------------------------------
        # edgeR
        # ----------------------------------------------
        with tempfile.TemporaryDirectory() as tmp:

            cf = os.path.join(tmp, "counts.csv")
            mf = os.path.join(tmp, "meta.csv")
            rf = os.path.join(tmp, "run.R")
            of = os.path.join(tmp, "out.csv")

            pd.DataFrame(
                counts.T,
                index=genes,
                columns=meta["sample"],
            ).to_csv(cf)

            meta.to_csv(mf, index=False)

            with open(rf, "w") as f:
                f.write(R_SCRIPT)

            try:
                subprocess.run(
                    ["Rscript", rf, cf, mf, of, group1, group2],
                    check=True,
                )
            except Exception:
                continue

            if not os.path.exists(of):
                continue

            res = pd.read_csv(of)

        # ----------------------------------------------
        # format like old code
        # ----------------------------------------------
        res = res.rename(columns={"PValue": "pval", "FDR": "fdr"})

        res["n_patients_group1"] = int(g1)
        res["n_patients_group2"] = int(g2)

        if isinstance(keys, tuple):
            for col, val in zip(split_by, keys):
                res[col] = val
            name = "_".join(map(str, keys))
        else:
            res[split_by[0]] = keys
            name = str(keys)

        if "fdr" in res.columns and "logFC" in res.columns:
            res = res.sort_values(
                ["fdr", "logFC"],
                ascending=[True, False],
            )

        res.to_csv(
            os.path.join(output_dir, f"{name}_DE.csv"),
            index=False,
        )