import numpy as np, pandas as pd, scipy.sparse as sp, joblib
from sklearn.neighbors import NearestNeighbors
from collections import Counter
from .io import  save_adata

import joblib
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from .io import load_adata

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import (
    accuracy_score,
    classification_report,
    confusion_matrix,
)
from sklearn.linear_model import SGDClassifier


def annotate_cells(
    leiden_path,
    output_path,
    marker_dict,
    clip_percentiles=(1, 99),
):
    """
        Annote automatiquement les cellules à partir de gènes marqueurs.

        La fonction charge un objet AnnData clusterisé, calcule pour chaque type
        cellulaire un score moyen basé sur l'expression normalisée des gènes
        marqueurs fournis, attribue à chaque cellule le type ayant le score le
        plus élevé, puis sauvegarde le résultat.

        :param leiden_path: Chemin vers le fichier .h5ad contenant les données
            clusterisées.
        :type leiden_path: str | Path

        :param output_path: Chemin du fichier .h5ad annoté en sortie.
        :type output_path: str | Path

        :param marker_dict: Dictionnaire associant chaque type cellulaire à une
            liste de gènes marqueurs.
        :type marker_dict: dict[str, list[str]]

        :param clip_percentiles: Percentiles utilisés pour tronquer les valeurs
            extrêmes d'expression avant standardisation.
        :type clip_percentiles: tuple[int | float, int | float]

        :return: Objet AnnData avec la colonne `cell_type` ajoutée dans `obs`.
        :rtype: anndata.AnnData
        """
    adata = load_adata(leiden_path)
    adata.obs_names_make_unique()

    X = adata.X.tocsr() if sp.issparse(adata.X) else np.asarray(adata.X)
    n_cells = adata.n_obs

    scores = {}

    for celltype, genes in marker_dict.items():
        genes = [g for g in genes if g in adata.var_names]

        if len(genes) == 0:
            continue

        vals = []

        for gene in genes:
            gidx = adata.var_names.get_loc(gene)

            col = X[:, gidx]
            col = col.toarray().ravel() if sp.issparse(col) else np.asarray(col).ravel()
            col = col.astype(np.float32, copy=False)

            lo = np.percentile(col, clip_percentiles[0])
            hi = np.percentile(col, clip_percentiles[1])
            col = np.clip(col, lo, hi)

            std = col.std()
            if std < 1e-6:
                continue

            col = (col - col.mean()) / std
            vals.append(col)

        if len(vals) == 0:
            continue

        scores[celltype] = np.mean(vals, axis=0).astype(np.float32)

    if len(scores) == 0:
        adata.obs["cell_type"] = "Unknown"
        save_adata(adata, output_path)
        return adata

    score_df = pd.DataFrame(scores, index=adata.obs_names)

    adata.obs["cell_type"] = score_df.idxmax(axis=1).astype(str).values

    save_adata(adata, output_path)
    return adata

def filter_high_confidence(
    annotated_path,
    output_path,
    latent_key="X_scVI",
    label_key="cell_type",
    min_consistency=0.7,
    n_neighbors=15,
    batch_size=5000,
):
    """
    Filtre les cellules annotées avec une forte confiance locale.

    La fonction mesure, pour chaque cellule, la cohérence de son annotation
    avec celle de ses plus proches voisins dans l'espace latent. Seules les
    cellules non "Unknown" dont la cohérence dépasse le seuil défini sont
    conservées puis sauvegardées.

    :param annotated_path: Chemin vers le fichier .h5ad contenant les cellules
        annotées.
    :type annotated_path: str | Path

    :param output_path: Chemin du fichier .h5ad filtré en sortie.
    :type output_path: str | Path

    :param latent_key: Clé de `adata.obsm` contenant la représentation latente
        utilisée pour la recherche de voisins.
    :type latent_key: str

    :param label_key: Colonne de `adata.obs` contenant les annotations
        cellulaires.
    :type label_key: str

    :param min_consistency: Seuil minimal de cohérence locale requis pour
        conserver une cellule.
    :type min_consistency: float

    :param n_neighbors: Nombre de voisins utilisés pour estimer la cohérence.
    :type n_neighbors: int

    :param batch_size: Taille des batchs utilisée pour calculer les voisins
        par blocs.
    :type batch_size: int

    :return: Sous-ensemble AnnData contenant uniquement les cellules
        high-confidence.
    :rtype: anndata.AnnData
    """
    adata = load_adata(annotated_path)
    adata.obs_names_make_unique()

    X = adata.obsm[latent_key].astype(np.float32, copy=False)
    y = adata.obs[label_key].astype(str).values

    nn = NearestNeighbors(
        n_neighbors=n_neighbors,
        algorithm="auto",
        n_jobs=-1
    )
    nn.fit(X)

    consistency = np.zeros(len(y), dtype=np.float32)

    for start in range(0, len(y), batch_size):
        end = min(start + batch_size, len(y))

        neigh = nn.kneighbors(
            X[start:end],
            return_distance=False
        )

        for j, row in enumerate(neigh):
            vals = y[row]
            vals = vals[vals != "Unknown"]

            if len(vals) == 0:
                continue

            consistency[start + j] = (
                Counter(vals).most_common(1)[0][1] / len(vals)
            )

        print(f"processed {end}/{len(y)}")

    mask = (
        (adata.obs[label_key].astype(str) != "Unknown").values &
        (consistency >= min_consistency)
    )

    adata.obs["consistency"] = consistency
    adata.obs["high_confidence"] = mask

    out = adata[mask]
    save_adata(out, output_path)

    return out


def train_classifier(
    highconf_path,
    model_path,
    latent_key="X_scVI",
    label_key="cell_type",
    test_size=0.2,
    random_state=0,
):
    """
    Version optimisée RAM :
    - évite les copies inutiles
    - libère adata tôt
    - conserve la logique initiale
    """

    # =====================================================
    # LOAD
    # =====================================================
    adata = load_adata(highconf_path)

    # copy=False => pas de copie si déjà float32
    X = adata.obsm[latent_key].astype(np.float32, copy=False)

    # extraction minimale labels
    y = adata.obs[label_key].to_numpy(dtype=str, copy=False)

    print("Data shape:", X.shape)

    # on n'a plus besoin du conteneur AnnData
    del adata

    # =====================================================
    # LABEL ENCODING
    # =====================================================
    le = LabelEncoder()
    y_enc = le.fit_transform(y)

    # libère y texte
    del y

    # =====================================================
    # SPLIT
    # =====================================================
    Xtr, Xte, ytr, yte = train_test_split(
        X,
        y_enc,
        test_size=test_size,
        stratify=y_enc,
        random_state=random_state
    )

    # libère matrices complètes
    del X, y_enc

    print("Train:", Xtr.shape, "| Test:", Xte.shape)

    # =====================================================
    # MODEL
    # =====================================================
    clf = SGDClassifier(
        loss="log_loss",
        penalty="l2",
        alpha=1e-4,
        max_iter=1000,
        tol=1e-3,
        n_jobs=-1,
        random_state=random_state,
        verbose=1
    )

    print("Training SGD logistic classifier...")
    clf.fit(Xtr, ytr)

    # =====================================================
    # EVAL
    # =====================================================
    pred = clf.predict(Xte)

    acc = accuracy_score(yte, pred)
    print("\nAccuracy:", acc)
    print("\nClassification report:\n")
    print(classification_report(yte, pred, target_names=le.classes_))

    # =====================================================
    # CONFUSION MATRIX EXPORT
    # =====================================================
    cm = confusion_matrix(yte, pred)

    fig_w = max(8, len(le.classes_) * 0.8)
    fig_h = max(6, len(le.classes_) * 0.6)

    plt.figure(figsize=(fig_w, fig_h))
    plt.imshow(cm, aspect="auto")
    plt.colorbar()

    plt.xticks(np.arange(len(le.classes_)), le.classes_, rotation=90)
    plt.yticks(np.arange(len(le.classes_)), le.classes_)

    plt.xlabel("Predicted")
    plt.ylabel("True")
    plt.title("Confusion Matrix")
    plt.tight_layout()

    out_dir = Path(model_path).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    cm_path = out_dir / "confusion_matrix.png"
    plt.savefig(cm_path, dpi=300, bbox_inches="tight")
    plt.close()

    print("Confusion matrix saved:", cm_path)

    # =====================================================
    # SAVE MODEL
    # =====================================================
    joblib.dump(
        {
            "model": clf,
            "label_encoder": le
        },
        model_path
    )

    print("Saved model:", model_path)

    return clf, acc

def predict_all_cells(
    latent_path,
    model_path,
    output_path,
    latent_key="X_scVI",
    pred_key="cell_type_model",
    proba_key="cell_type_model_score",
):
    """
    Prédit les annotations cellulaires sur l'ensemble des cellules à partir
    d'un classifieur entraîné.

    La fonction charge un objet AnnData contenant la représentation latente,
    charge le modèle supervisé sauvegardé, prédit les labels cellulaires pour
    chaque cellule, calcule les probabilités associées lorsque disponibles,
    ajoute les résultats dans `adata.obs`, puis sauvegarde l'objet annoté.

    :param latent_path: Chemin vers le fichier .h5ad contenant les données
        et la représentation latente.
    :type latent_path: str | Path

    :param model_path: Chemin vers le fichier du modèle supervisé sauvegardé.
    :type model_path: str | Path

    :param output_path: Chemin du fichier .h5ad annoté en sortie.
    :type output_path: str | Path

    :param latent_key: Clé de `adata.obsm` contenant l'embedding utilisé
        pour la prédiction.
    :type latent_key: str

    :param pred_key: Nom de la colonne ajoutée dans `adata.obs` contenant
        les labels prédits.
    :type pred_key: str

    :param proba_key: Nom de la colonne ajoutée dans `adata.obs` contenant
        le score de confiance ou la probabilité maximale prédite.
    :type proba_key: str

    :return: Objet AnnData enrichi avec les prédictions.
    :rtype: anndata.AnnData
    """


    # =====================================================
    # LOAD DATA
    # =====================================================
    adata = load_adata(latent_path)
    X = adata.obsm[latent_key].astype(np.float32, copy=False)

    # =====================================================
    # LOAD MODEL
    # =====================================================
    obj = joblib.load(model_path)
    clf = obj["model"]
    le = obj["label_encoder"]

    # =====================================================
    # PREDICT
    # =====================================================
    pred_enc = clf.predict(X)
    pred = le.inverse_transform(pred_enc)

    adata.obs[pred_key] = pred

    # =====================================================
    # CONFIDENCE SCORE
    # =====================================================
    if hasattr(clf, "predict_proba"):
        probs = clf.predict_proba(X)
        score = probs.max(axis=1)

    elif hasattr(clf, "decision_function"):
        dec = clf.decision_function(X)

        if dec.ndim == 1:
            score = 1 / (1 + np.exp(-dec))
        else:
            exp = np.exp(dec - dec.max(axis=1, keepdims=True))
            soft = exp / exp.sum(axis=1, keepdims=True)
            score = soft.max(axis=1)
    else:
        score = np.ones(adata.n_obs)

    adata.obs[proba_key] = score.astype(np.float32)

    # =====================================================
    # SAVE
    # =====================================================
    save_adata(adata, output_path)

    print("Annotated cells:", adata.n_obs)
    print("Saved:", output_path)

    return adata