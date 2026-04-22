import scanpy as sc
from sklearn.cluster import MiniBatchKMeans
import numpy as np

from .io import load_adata, save_adata
from .utils import ensure_dir

def compute_leiden(latent_path, output_path, resolution=1.0, key_added='leiden'):
    """
    Effectue un clustering de Leiden à partir d'une représentation latente.

    La fonction charge un objet AnnData contenant un embedding (ex. `X_scVI`),
    calcule le graphe de voisinage, applique l'algorithme de Leiden puis
    sauvegarde le résultat annoté.

    :param latent_path: Chemin vers le fichier .h5ad contenant les données
        et la représentation latente.
    :type latent_path: str | Path

    :param output_path: Chemin du fichier de sortie .h5ad.
    :type output_path: str | Path

    :param resolution: Paramètre de résolution du clustering Leiden.
        Plus la valeur est élevée, plus le nombre de clusters tend à augmenter.
    :type resolution: float

    :param key_added: Nom de la colonne ajoutée dans `adata.obs`
        contenant les labels de clusters.
    :type key_added: str

    :return: Objet AnnData avec les clusters calculés.
    :rtype: anndata.AnnData
    """
    adata=load_adata(latent_path)
    sc.pp.neighbors(adata, use_rep='X_scVI')
    sc.tl.leiden(adata, resolution=resolution, key_added=key_added)
    save_adata(adata, output_path)
    return adata

def subsample_umap(
    latent_path,
    output_png,
    group_key,
    groups=None,
    n_cells=8000,
    random_state=0,
    basis="X_scVI",
):
    """
    Génère une visualisation UMAP à partir d'un sous-échantillonnage équilibré
    des cellules.

    La fonction charge un objet AnnData, filtre éventuellement certaines
    catégories, sélectionne un nombre équilibré de cellules par groupe,
    calcule un UMAP sur ce sous-ensemble puis sauvegarde la figure au format PNG.

    :param latent_path: Chemin vers le fichier .h5ad contenant les données
        et la représentation latente.
    :type latent_path: str | Path

    :param output_png: Chemin du fichier image PNG de sortie.
    :type output_png: str | Path

    :param group_key: Colonne de `adata.obs` utilisée pour colorer les cellules
        sur l'UMAP.
    :type group_key: str

    :param groups: Valeurs spécifiques de `group_key` à conserver.
        Si None, toutes les catégories sont utilisées.
    :type groups: list[str] | str | None

    :param n_cells: Nombre maximal de cellules à afficher après
        sous-échantillonnage.
    :type n_cells: int

    :param random_state: Graine aléatoire pour la reproductibilité du
        sous-échantillonnage.
    :type random_state: int

    :param basis: Clé de `adata.obsm` contenant l'embedding utilisé
        pour calculer les voisins (ex. `X_scVI`).
    :type basis: str

    :return: None
    :rtype: None
    """
    import numpy as np
    import scanpy as sc
    import matplotlib.pyplot as plt
    from pathlib import Path

    # --------------------------------------------------
    # Load
    # --------------------------------------------------
    adata = load_adata(latent_path)

    if group_key not in adata.obs.columns:
        raise ValueError(f"{group_key} absent de adata.obs")

    if basis not in adata.obsm:
        raise ValueError(f"{basis} absent de adata.obsm")

    obs_values = adata.obs[group_key].astype(str)

    # --------------------------------------------------
    # Optional filtering
    # --------------------------------------------------
    if groups is not None:
        if isinstance(groups, str):
            groups = [groups]

        mask = obs_values.isin(groups)
        adata = adata[mask].copy()
        obs_values = adata.obs[group_key].astype(str)

    if adata.n_obs == 0:
        raise ValueError("Aucune cellule après filtrage")

    # --------------------------------------------------
    # Balanced subsampling
    # --------------------------------------------------
    rng = np.random.default_rng(random_state)

    cats = np.array(obs_values.unique())
    per_group = max(1, n_cells // len(cats))

    idx = []
    values = obs_values.values

    for g in cats:
        g_idx = np.where(values == g)[0]
        take = min(len(g_idx), per_group)
        idx.extend(rng.choice(g_idx, size=take, replace=False))

    idx = np.array(idx)

    sub = adata[idx].copy()
    del adata

    # --------------------------------------------------
    # Compute UMAP on subset only
    # --------------------------------------------------
    sc.pp.neighbors(
        sub,
        use_rep=basis,
        n_neighbors=10,
    )

    sc.tl.umap(sub)

    # --------------------------------------------------
    # Save figure
    # --------------------------------------------------
    output_png = Path(output_png)
    output_png.parent.mkdir(parents=True, exist_ok=True)

    sc.pl.umap(
        sub,
        color=group_key,
        show=False,
        frameon=False,
        size=8,
    )

    plt.savefig(output_png, dpi=250, bbox_inches="tight")
    plt.close()

    del sub