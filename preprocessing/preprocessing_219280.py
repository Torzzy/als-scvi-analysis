import os
import scanpy as sc
import numpy as np
from pathlib import Path


def preprocess_h5_scvi(
        input_path,
        output_dir,
        sample_name=None,
        condition=None,
        region=None,
        patient_id=None,
):
    """
    Prétraite un fichier 10x Genomics `.h5` pour une utilisation avec scVI.

    La fonction lit une matrice d'expression au format 10x, rend les noms de
    gènes uniques, ajoute les métadonnées d'échantillon, calcule des métriques
    QC minimales, prépare la couche de counts brute pour scVI, convertit les
    données au bon format, puis sauvegarde un fichier `.h5ad`.

    :param input_path: Chemin vers le fichier d'entrée `.h5`.
    :type input_path: str | Path

    :param output_dir: Dossier de sortie où sauvegarder le fichier `.h5ad`.
    :type output_dir: str | Path

    :param sample_name: Nom de l'échantillon à enregistrer dans `obs`.
        Si None, utilise le nom du fichier.
    :type sample_name: str | None

    :param condition: Condition biologique associée à l'échantillon
        (ex. ALS, CTRL).
    :type condition: str | None

    :param region: Région anatomique associée à l'échantillon.
    :type region: str | None

    :param patient_id: Identifiant patient associé à l'échantillon.
    :type patient_id: str | None

    :return: None
    :rtype: None
    """
    print(f" Lecture : {input_path}")

    adata = sc.read_10x_h5(input_path)
    adata.var_names_make_unique()

    # =========================
    # METADATA
    # =========================
    adata.obs['sample'] = sample_name if sample_name else os.path.basename(input_path)
    adata.obs['condition'] = condition if condition else "unknown"
    adata.obs['region'] = region if region else "unknown"
    adata.obs['patient_id'] = patient_id if patient_id else "unknown"

    dataset_name = Path(input_path).parent.name
    adata.obs['batch'] = dataset_name
    adata.obs['dataset'] = dataset_name

    # =========================
    # MINIMAL METRICS
    # =========================
    X = adata.X

    adata.obs['n_counts'] = np.asarray(X.sum(axis=1)).ravel()

    if hasattr(X, "getnnz"):
        adata.obs['n_genes'] = X.getnnz(axis=1)
    else:
        adata.obs['n_genes'] = np.count_nonzero(X, axis=1)

    # mitochondrial genes flag
    adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')

    # =========================
    #  scVI
    # =========================
    adata.layers["counts"] = adata.X.copy()

    if not isinstance(adata.X, np.ndarray):
        adata.X = adata.X.tocsr()

    adata.X = adata.X.astype('float32')

    # =========================
    # SAVE
    # =========================
    os.makedirs(output_dir, exist_ok=True)

    output_name = os.path.basename(input_path).replace(".h5", ".h5ad")
    output_path = os.path.join(output_dir, output_name)

    adata.write(output_path, compression="gzip")

    print(" Terminé\n")

    del adata


def batch_preprocess_scvi_auto(
    input_dir,
    output_root
):
    """
    Prétraite automatiquement tous les fichiers `.h5` d'un dossier.

    La fonction parcourt un répertoire d'entrée, détecte les fichiers 10x,
    infère automatiquement la condition biologique, la région anatomique,
    l'identifiant patient et le nom d'échantillon à partir du nom des fichiers,
    exclut les échantillons FTD, puis lance le prétraitement individuel de
    chaque dataset.

    Les résultats sont sauvegardés dans :
    `output_root / nom_du_dossier_input /`

    :param input_dir: Dossier contenant les fichiers `.h5` à traiter.
    :type input_dir: str | Path

    :param output_root: Dossier racine de sortie.
    :type output_root: str | Path

    :return: None
    :rtype: None
    """

    input_dir = Path(input_dir)
    output_root = Path(output_root)

    dataset_name = input_dir.name
    output_dir = output_root / dataset_name
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f" Input : {input_dir}")
    print(f" Output : {output_dir}\n")

    files = list(input_dir.glob("*.h5"))
    print(f" {len(files)} fichiers détectés\n")

    for f in files:
        fname = f.name

        # =========================
        # EXCLUSION FTD
        # =========================
        if "FTD" in fname:
            print(f" Skip (FTD) : {fname}")
            continue

        # =========================
        # CONDITION
        # =========================
        if "ALS" in fname:
            condition = "ALS"
        elif "Control" in fname:
            condition = "CTRL"
        else:
            print(f" Condition inconnue : {fname}")
            condition = "unknown"

        # =========================
        # REGION
        # =========================
        if "MotorCortex" in fname:
            region = "MCX"
        elif "FrontalCortex" in fname:
            region = "FCX"
        else:
            print(f" Région inconnue : {fname}")
            region = "unknown"

        # =========================
        # SAMPLE NAME
        # =========================
        sample_name = fname.replace(".h5", "")

        # =========================
        # PATIENT ID
        # =========================
        patient_id = sample_name.split("_")[0]

        print(f"\n Processing : {fname}")
        print(f"   condition = {condition}, region = {region}")

        preprocess_h5_scvi(
            input_path=str(f),
            output_dir=str(output_dir),
            sample_name=sample_name,
            condition=condition,
            region=region,
            patient_id=patient_id
        )

    print(f"\n Terminé : résultats dans {output_dir}")