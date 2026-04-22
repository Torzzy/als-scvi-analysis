import numpy as np, scvi
from .io import load_adata, save_adata

def compute_scvi_latent(h5ad_path, model_path, output_path, batch_size=1024):
    """
    Calcule la représentation latente scVI d'un jeu de données single-cell.

    La fonction charge un objet AnnData ainsi qu'un modèle scVI entraîné,
    projette les cellules dans l'espace latent, stocke cet embedding dans
    `adata.obsm["X_scVI"]`, puis sauvegarde le résultat.

    :param h5ad_path: Chemin vers le fichier .h5ad contenant les données brutes.
    :type h5ad_path: str | Path

    :param model_path: Chemin vers le modèle scVI sauvegardé.
    :type model_path: str | Path

    :param output_path: Chemin du fichier .h5ad de sortie contenant
        la représentation latente.
    :type output_path: str | Path

    :param batch_size: Taille des batchs utilisée pour calculer
        l'embedding latent.
    :type batch_size: int

    :return: Objet AnnData enrichi avec `X_scVI`.
    :rtype: anndata.AnnData
    """
    adata=load_adata(h5ad_path)
    model=scvi.model.SCVI.load(model_path, adata=adata)
    latent=model.get_latent_representation(adata=adata,batch_size=batch_size).astype(np.float32)
    adata.obsm['X_scVI']=latent
    save_adata(adata, output_path)
    return adata