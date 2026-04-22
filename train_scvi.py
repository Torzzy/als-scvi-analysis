import scanpy as sc
import scvi
import numpy as np
import os
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report

from utils import get_next_run_dir

def train_scvi(
    h5ad_path,
    batch_key="patient_id",
    n_epochs=100,
    batch_size=1024,
):
    """
    Entraîne un modèle scVI à partir d'un fichier AnnData (.h5ad).

    Charge les données single-cell, vérifie les identifiants de cellules,
    configure scVI avec les covariables disponibles, entraîne le modèle,
    puis sauvegarde les poids dans un nouveau dossier de run.

    :param h5ad_path: Chemin vers le fichier .h5ad contenant les données d'entrée.
    :type h5ad_path: str | Path

    :param batch_key: Nom de la colonne `adata.obs` utilisée comme batch biologique
        ou technique (ex: patient_id).
    :type batch_key: str

    :param n_epochs: Nombre maximal d'époques d'entraînement.
    :type n_epochs: int

    :param batch_size: Taille des mini-batchs utilisée pendant l'entraînement.
    :type batch_size: int

    :return: Modèle scVI entraîné.
    :rtype: scvi.model.SCVI
    """

    output_dir = get_next_run_dir()

    print("📥 Loading dataset...")
    adata = sc.read_h5ad(h5ad_path)

    print("📊 Shape:", adata.shape)

    # -------------------------
    # 🔧 FIX noms cellules
    # -------------------------
    if not adata.obs_names.is_unique:
        print("⚠️ Fixing duplicate cell names...")
        adata.obs_names_make_unique()

    # -------------------------
    # 🔍 CHECK batch
    # -------------------------
    if batch_key not in adata.obs.columns:
        raise ValueError(f"❌ {batch_key} not found in obs")

    print(f"📊 Batch ({batch_key}) unique values:",
          adata.obs[batch_key].nunique())

    # -------------------------
    # ⚙️ SETUP SCVI
    # -------------------------
    print("⚙️ Setting up scVI...")

    scvi.model.SCVI.setup_anndata(
        adata,
        batch_key=batch_key,
        categorical_covariate_keys=[
            col for col in ["condition", "region"]
            if col in adata.obs.columns
        ]
    )

    # -------------------------
    # 🚀 MODEL
    # -------------------------
    print("🚀 Initializing model...")

    model = scvi.model.SCVI(
        adata,
        n_layers=2,
        n_latent=30,
        gene_likelihood="nb"
    )

    # -------------------------
    # 🧪 TRAIN
    # -------------------------
    print("🧪 Training...")

    model.train(
        max_epochs=n_epochs,
        batch_size=batch_size
    )

    # -------------------------
    # 💾 SAVE MODEL
    # -------------------------
    print("💾 Saving model...")

    model_dir = os.path.join(output_dir, "model")
    model.save(model_dir, overwrite=True)


    return model

# TEST ENTRAINEMENT MODEL
def evaluate_region_classification(
    adata,
    latent_key="X_scVI",
    label_key="region",
    test_size=0.2,
    random_state=42,
    max_iter=1000
):
    """
    Évalue la capacité de l'espace latent à prédire la région d'origine
    des cellules via une régression logistique supervisée.

    Les embeddings contenus dans `adata.obsm[latent_key]` sont utilisés
    comme variables explicatives, tandis que `adata.obs[label_key]`
    fournit les labels de classes. Les données sont séparées en jeu
    d'entraînement et de test avant apprentissage puis évaluation.

    :param adata: Objet AnnData contenant embeddings et annotations.
    :type adata: anndata.AnnData

    :param latent_key: Clé de `adata.obsm` contenant la représentation latente.
    :type latent_key: str

    :param label_key: Colonne de `adata.obs` contenant les labels à prédire.
    :type label_key: str

    :param test_size: Proportion des données utilisée pour le jeu de test.
    :type test_size: float

    :param random_state: Graine aléatoire pour la reproductibilité du split.
    :type random_state: int

    :param max_iter: Nombre maximal d'itérations de la régression logistique.
    :type max_iter: int

    :return: Tuple contenant l'accuracy sur le jeu de test et le rapport
        de classification complet.
    :rtype: tuple[float, str]
    """

    # =========================
    #  Extract data
    # =========================
    X = adata.obsm[latent_key]
    y = adata.obs[label_key].values

    # remove NaNs if any
    mask = ~np.isnan(X).any(axis=1)
    X = X[mask]
    y = y[mask]

    # =========================
    #  Train / Test split
    # =========================
    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=test_size,
        random_state=random_state,
        stratify=y
    )

    # =========================
    #  Model
    # =========================
    clf = LogisticRegression(
        max_iter=max_iter,
        n_jobs=-1
    )

    clf.fit(X_train, y_train)

    # =========================
    #  Evaluation
    # =========================
    y_pred = clf.predict(X_test)

    acc = accuracy_score(y_test, y_pred)
    report = classification_report(y_test, y_pred)

    # =========================
    #  Print clean
    # =========================
    print("\n📊 Logistic Regression (region prediction)")
    print(f"Accuracy: {acc:.4f}\n")
    print(report)

    return acc, report

def compute_umap_from_scvi(
    h5ad_path,
    model_path,
    output_dir,
    batch_size=1024,
    n_cells_umap=100000
):
    """
    Calcule un embedding latent scVI puis génère un UMAP sur un sous-échantillon
    de cellules.

    La fonction charge les données et un modèle scVI préalablement entraîné,
    projette les cellules dans l'espace latent `X_scVI`, sous-échantillonne
    éventuellement les observations pour accélérer l'UMAP, calcule le graphe
    de voisinage puis sauvegarde l'objet AnnData réduit.

    :param h5ad_path: Chemin vers le fichier .h5ad d'entrée.
    :type h5ad_path: str | Path

    :param model_path: Chemin vers le dossier du modèle scVI sauvegardé.
    :type model_path: str | Path

    :param output_dir: Dossier dans lequel sauvegarder le résultat UMAP.
    :type output_dir: str | Path

    :param batch_size: Taille des batchs utilisée pour calculer
        la représentation latente.
    :type batch_size: int

    :param n_cells_umap: Nombre maximal de cellules conservées pour
        le calcul de l'UMAP.
    :type n_cells_umap: int

    :return: None
    :rtype: None
    """
    print(" Loading dataset...")
    adata = sc.read_h5ad(h5ad_path)

    print(" Loading model...")
    model = scvi.model.SCVI.load(model_path, adata=adata)

    # -------------------------
    #  LATENT (BATCHED)
    # -------------------------
    print(" Computing latent (batched)...")

    latent = model.get_latent_representation(batch_size=batch_size)
    latent = latent.astype(np.float32)

    adata.obsm["X_scVI"] = latent
    del latent

    # -------------------------
    #  SUBSAMPLE
    # -------------------------
    print(f" Subsampling to {n_cells_umap} cells for UMAP...")

    if adata.n_obs > n_cells_umap:
        idx = np.random.choice(adata.n_obs, n_cells_umap, replace=False)
        adata_sub = adata[idx].copy()
    else:
        adata_sub = adata.copy()

    # -------------------------
    #  NEIGHBORS
    # -------------------------
    print(" Computing neighbors...")

    sc.pp.neighbors(
        adata_sub,
        use_rep="X_scVI",
        n_neighbors=10,
        method="umap"
    )

    # -------------------------
    #  UMAP
    # -------------------------
    print(" Computing UMAP...")

    sc.tl.umap(adata_sub)

    # -------------------------
    #  PLOT
    # -------------------------
    print(" Plotting...")

    for col in ["batch", "dataset", "patient_id", "condition"]:
        if col in adata_sub.obs.columns:
            sc.pl.umap(adata_sub, color=col, show=True)

    # -------------------------
    #  SAVE FIG
    # -------------------------
    print(" Saving UMAP object...")

    out_path = f"{output_dir}/umap_subset.h5ad"
    adata_sub.write(out_path, compression="gzip")

    print(" DONE")