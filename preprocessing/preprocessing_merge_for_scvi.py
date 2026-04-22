
import numpy as np
import scanpy as sc
from pathlib import Path
import anndata as ad
import tempfile
import shutil


def harmonize_gene_names(adata):
    """
        Harmonise les noms de gènes d'un objet AnnData.

        La fonction supprime les suffixes de version des identifiants Ensembl,
        remplace les identifiants ENSG par les symboles géniques lorsque la
        colonne `common_name` est disponible, standardise l'écriture des noms,
        retire les gènes vides, corrige certains problèmes d'index AnnData,
        puis rend les noms de variables uniques.

        :param adata: Objet AnnData contenant les gènes dans `var_names`.
        :type adata: anndata.AnnData

        :return: Objet AnnData avec noms de gènes harmonisés.
        :rtype: anndata.AnnData
        """
    # enlever version ENSG
    adata.var_names = adata.var_names.str.replace(r"\.\d+$", "", regex=True)

    # ENSG → SYMBOL si possible
    if adata.var_names.str.startswith("ENSG").mean() > 0.5:
        if "common_name" in adata.var.columns:
            adata.var_names = adata.var["common_name"].astype(str)

    # nettoyage
    adata.var_names = adata.var_names.str.upper().str.strip()

    # enlever gènes vides
    mask = adata.var_names != ""
    adata = adata[:, mask]

    # fix bug anndata
    if adata.var.index.name == "common_name":
        adata.var.index.name = None

    if "common_name" in adata.var.columns:
        adata.var = adata.var.drop(columns=["common_name"])

    # rendre unique (CRUCIAL)
    adata.var_names_make_unique()

    return adata


def build_common_gene_datasets(
    dataset_dirs,
    output_root="datasets_common_genes"
):
    """
    Construit des datasets harmonisés partageant le même ensemble de gènes.

    La fonction parcourt plusieurs dossiers contenant des fichiers `.h5ad`,
    harmonise les noms de gènes de chaque fichier, calcule l'intersection
    globale des gènes communs, puis sauvegarde chaque dataset restreint à cet
    ensemble commun tout en conservant la structure des dossiers d'origine.

    Le traitement est réalisé fichier par fichier afin de limiter l'usage
    mémoire.

    :param dataset_dirs: Liste des dossiers contenant les fichiers `.h5ad`
        à harmoniser.
    :type dataset_dirs: list[str | Path]

    :param output_root: Dossier racine de sortie contenant les datasets
        filtrés sur les gènes communs.
    :type output_root: str | Path

    :return: None
    :rtype: None
    """

    print(" Collecting files...")

    dataset_dirs = [Path(d) for d in dataset_dirs]

    all_files = []
    for d in dataset_dirs:
        files = list(d.glob("*.h5ad"))
        all_files.extend(files)

    print(f" Total files: {len(all_files)}")

    # -------------------------
    #  INTERSECTION DES GÈNES
    # -------------------------
    print(" Computing gene intersection...")

    gene_sets = []

    for i, f in enumerate(all_files):
        print(f"   [{i+1}/{len(all_files)}] {f.name}")

        adata = sc.read_h5ad(f)
        adata = harmonize_gene_names(adata)

        gene_sets.append(set(adata.var_names))

        del adata  #  libère RAM

    common_genes = sorted(set.intersection(*gene_sets))

    print(f"\n Common genes: {len(common_genes)}")

    if len(common_genes) == 0:
        raise ValueError(" ZERO gene overlap → mapping encore cassé")

    # -------------------------
    #  SAVE DATASETS
    # -------------------------
    print("\n Saving filtered datasets...")

    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    for d in dataset_dirs:

        dataset_name = d.name
        out_dir = output_root / dataset_name
        out_dir.mkdir(parents=True, exist_ok=True)

        files = list(d.glob("*.h5ad"))

        print(f"\n Dataset: {dataset_name} ({len(files)} files)")

        for i, f in enumerate(files):
            print(f"   [{i+1}/{len(files)}] {f.name}")

            adata = sc.read_h5ad(f)
            adata = harmonize_gene_names(adata)

            #  SUBSET STRICT (PAS de copy inutile)
            adata = adata[:, common_genes]

            # sécurité noms
            adata.obs_names_make_unique()
            adata.var_names_make_unique()

            # save
            out_path = out_dir / f.name
            adata.write(out_path, compression="gzip")

            del adata  #  libère RAM

    print("\n DONE")
    print(f" Output: {output_root}")

def qc_hvg_global_from_root(
    input_root,
    output_root="datasets_scvi_ready",
    min_genes=200,
    max_mt=0.2,
    min_cells_per_gene=3,
    n_hvg=3000
):
    """
    Applique un filtrage QC global et une sélection HVG sur plusieurs datasets.

    La fonction parcourt tous les fichiers `.h5ad` contenus dans un dossier
    racine et exécute trois passes successives :

    1. Calcule les statistiques globales d'expression pour sélectionner les
       gènes hautement variables (HVG).
    2. Applique un contrôle qualité cellulaire global et compte les gènes
       détectés après filtrage.
    3. Sauvegarde les datasets finaux filtrés avec une structure identique
       à celle des dossiers d'entrée.

    Les fichiers générés sont prêts pour une utilisation avec scVI.

    :param input_root: Dossier racine contenant des sous-dossiers avec les
        fichiers `.h5ad`.
    :type input_root: str | Path

    :param output_root: Dossier racine de sortie.
    :type output_root: str | Path

    :param min_genes: Nombre minimal de gènes détectés par cellule pour
        conserver une cellule.
    :type min_genes: int

    :param max_mt: Fraction maximale d'expression mitochondriale autorisée
        par cellule.
    :type max_mt: float

    :param min_cells_per_gene: Nombre minimal de cellules exprimant un gène
        pour conserver ce gène.
    :type min_cells_per_gene: int

    :param n_hvg: Nombre de gènes hautement variables à sélectionner.
    :type n_hvg: int

    :return: None
    :rtype: None
    """

    input_root = Path(input_root)
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    all_files = list(input_root.glob("*/*.h5ad"))

    print(f" Total files: {len(all_files)}")

    # =========================================================
    #  PASS 1 — stats globales pour HVG (STREAMING)
    # =========================================================
    print(" Pass 1: computing global gene stats...")

    gene_sum = None
    gene_sq_sum = None
    total_cells = 0

    for i, f in enumerate(all_files):
        print(f"  ➤ [{i+1}/{len(all_files)}] {f.name}")

        adata = sc.read_h5ad(f)

        X = adata.X

        # somme
        s = np.asarray(X.sum(axis=0)).ravel()

        # somme des carrés
        if hasattr(X, "multiply"):
            sq = np.asarray(X.multiply(X).sum(axis=0)).ravel()
        else:
            sq = np.asarray((X**2).sum(axis=0)).ravel()

        if gene_sum is None:
            gene_sum = s
            gene_sq_sum = sq
        else:
            gene_sum += s
            gene_sq_sum += sq

        total_cells += adata.n_obs

        del adata

    mean = gene_sum / total_cells
    var = gene_sq_sum / total_cells - mean**2

    # HVG selection
    top_idx = np.argsort(var)[-n_hvg:]
    hvg_mask = np.zeros(len(var), dtype=bool)
    hvg_mask[top_idx] = True

    print(f" HVG selected: {hvg_mask.sum()}")

    # =========================================================
    #  PASS 2 — QC GLOBAL + comptage gènes
    # =========================================================
    print("\n Pass 2: QC + gene detection...")

    gene_detect_count = np.zeros(len(var))

    qc_info = []

    for i, f in enumerate(all_files):
        print(f"   [{i+1}/{len(all_files)}] {f.name}")

        adata = sc.read_h5ad(f)

        X = adata.X

        n_counts = np.asarray(X.sum(axis=1)).ravel()
        n_genes = X.getnnz(axis=1) if hasattr(X, "getnnz") else np.count_nonzero(X, axis=1)

        mt_mask = adata.var_names.str.startswith("MT-")

        if mt_mask.sum() > 0:
            mt_counts = np.asarray(adata[:, mt_mask].X.sum(axis=1)).ravel()
            pct_mt = mt_counts / n_counts
        else:
            pct_mt = np.zeros_like(n_counts)

        keep_cells = (
            (n_genes > min_genes) &
            (pct_mt < max_mt)
        )

        # comptage gènes (post-QC)
        if hasattr(X, "getnnz"):
            gene_detect_count += np.asarray(X[keep_cells].getnnz(axis=0)).ravel()
        else:
            gene_detect_count += np.count_nonzero(X[keep_cells], axis=0)

        qc_info.append((f, keep_cells))

        del adata

    # filtre gènes globaux
    gene_keep_mask = gene_detect_count >= min_cells_per_gene

    # combine HVG + gene filter
    final_gene_mask = hvg_mask & gene_keep_mask

    print(f"\n Genes kept after filter: {final_gene_mask.sum()}")

    # =========================================================
    #  PASS 3 — SAVE FINAL FILES
    # =========================================================
    print("\n Pass 3: saving final datasets...")

    for i, (f, keep_cells) in enumerate(qc_info):
        print(f"   [{i+1}/{len(qc_info)}] {f.name}")

        adata = sc.read_h5ad(f)

        # QC cells
        adata = adata[keep_cells]

        # genes
        adata = adata[:, final_gene_mask]

        # scVI
        adata.layers["counts"] = adata.X.copy()

        if not isinstance(adata.X, np.ndarray):
            adata.X = adata.X.tocsr()

        adata.X = adata.X.astype("float32")

        # save structure identique
        dataset_name = f.parent.name
        out_dir = output_root / dataset_name
        out_dir.mkdir(parents=True, exist_ok=True)

        out_path = out_dir / f.name
        adata.write(out_path, compression="gzip")

        del adata

    print("\n DONE")
    print(f" Output: {output_root}")



def merge_scvi_ready_datasets(
    input_root,
    output_path="final_scvi.h5ad"
):
    """
    Fusionne plusieurs datasets `.h5ad` préparés pour scVI en un seul fichier.

    La fonction recherche tous les fichiers `.h5ad` dans les sous-dossiers du
    dossier racine, applique des vérifications minimales, crée des fichiers
    temporaires harmonisés, puis effectue une concaténation sur disque avec
    intersection des gènes communs (`join="inner"`). Les fichiers temporaires
    sont ensuite supprimés.

    Cette approche limite l'utilisation mémoire pour les gros volumes de
    données.

    :param input_root: Dossier racine contenant les sous-dossiers avec les
        fichiers `.h5ad` à fusionner.
    :type input_root: str | Path

    :param output_path: Chemin du fichier `.h5ad` fusionné en sortie.
    :type output_path: str | Path

    :return: Chemin du fichier fusionné généré.
    :rtype: str | Path
    """

    input_root = Path(input_root)
    files = list(input_root.glob("*/*.h5ad"))

    print(f" Total files: {len(files)}")

    if len(files) == 0:
        raise ValueError(" Aucun fichier trouvé")

    # -------------------------------------------------
    #  dossier temporaire isolé
    # -------------------------------------------------
    tmp_dir = Path(tempfile.mkdtemp(prefix="scvi_merge_"))
    print(f" Temp dir: {tmp_dir}")

    tmp_files = []

    # -------------------------------------------------
    #  préparation fichiers
    # -------------------------------------------------
    print(" Preparing files...")

    for i, f in enumerate(files):
        print(f"   [{i+1}/{len(files)}] {f.name}")

        adata = sc.read_h5ad(f)

        # sécurité minimale
        adata.obs_names_make_unique()
        adata.var_names_make_unique()

        if "batch" not in adata.obs:
            adata.obs["batch"] = f.stem

        tmp_f = tmp_dir / f.name
        adata.write(tmp_f)

        tmp_files.append(tmp_f)

        del adata

    # -------------------------------------------------
    #  merge disk
    # -------------------------------------------------
    print("\n Merging on disk...")

    ad.experimental.concat_on_disk(
        tmp_files,
        output_path,
        join="inner"
    )

    print(" Merge done")

    # -------------------------------------------------
    #  cleanup
    # -------------------------------------------------
    print(" Cleaning temp dir...")
    shutil.rmtree(tmp_dir)

    print(" DONE")

    return output_path