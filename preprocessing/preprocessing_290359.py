import os
import scanpy as sc
import numpy as np
from pathlib import Path

def preprocess_all_patients(input_file, output_root):
    """
    Prétraite un fichier `.h5ad` global en séparant les données par patient.

    La fonction charge un dataset AnnData en mode lazy/backed, identifie les
    patients présents, extrait successivement les cellules de chaque patient en
    mémoire, ajoute les métadonnées harmonisées, calcule des métriques QC
    minimales, prépare la couche de counts brute pour scVI, optimise le format
    mémoire, puis sauvegarde un fichier `.h5ad` par patient.

    Les résultats sont enregistrés dans :
    `output_root / preprocessed_data / 290359 /`

    :param input_file: Chemin vers le fichier `.h5ad` source contenant tous
        les patients.
    :type input_file: str | Path

    :param output_root: Dossier racine de sortie.
    :type output_root: str | Path

    :return: None
    :rtype: None
    """

    print(f" Lecture lazy : {input_file}")
    adata = sc.read_h5ad(input_file, backed='r')

    patients = adata.obs['Sample ID'].unique()
    print(f" Total patients : {len(patients)}\n")

    output_root = Path(os.path.join(output_root,"preprocessed_data/290359"))
    output_root.mkdir(parents=True, exist_ok=True)

    for patient_id in patients:
        print(f" Traitement patient : {patient_id}")

        # Chargement en mémoire uniquement du patient
        adata_patient = adata[adata.obs['Sample ID'] == patient_id].to_memory()

        # =========================
        # METADATA
        # =========================
        adata_patient.obs['sample'] = patient_id
        adata_patient.obs['patient_id'] = patient_id
        adata_patient.obs['region'] = "FX"

        # CONDITION (ALS vs CTRL)
        disease = adata_patient.obs['Disease?'].astype(str)
        adata_patient.obs['condition'] = disease.map({
            'Yes': 'ALS',
            'No': 'CTRL',
        }).fillna('unknown')

        # batch et dataset
        adata_patient.obs['batch'] = "290359"
        adata_patient.obs['dataset'] = "290359"

        # =========================
        # MINIMAL METRICS (léger, sans filtrage)
        # =========================
        X = adata_patient.X

        adata_patient.obs['n_counts'] = np.asarray(X.sum(axis=1)).ravel()

        if hasattr(X, "getnnz"):
            adata_patient.obs['n_genes'] = X.getnnz(axis=1)
        else:
            adata_patient.obs['n_genes'] = np.count_nonzero(X, axis=1)

        # mitochondrial genes flag (utile plus tard)
        adata_patient.var['mt'] = adata_patient.var_names.str.upper().str.startswith('MT-')

        # =========================
        # IMPORTANT scVI
        # =========================
        adata_patient.layers["counts"] = adata_patient.X.copy()

        # =========================
        # OPTIM RAM
        # =========================
        if not isinstance(adata_patient.X, np.ndarray):
            adata_patient.X = adata_patient.X.tocsr()

        adata_patient.X = adata_patient.X.astype('float32')

        # =========================
        # SAVE
        # =========================
        output_path = output_root / f"{patient_id}.h5ad"
        adata_patient.write(output_path, compression="gzip")

        print(f" Sauvegardé : {output_path}\n")

        # libérer mémoire
        del adata_patient

    print(f" Terminé : tous les patients sauvegardés dans {output_root}")