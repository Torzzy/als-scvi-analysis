from pathlib import Path

from .config import DEFAULT_ANALYSIS_ROOT, FILES
from .utils import ensure_dir, reset_folder
from .latent import compute_scvi_latent
from .clustering import compute_leiden, subsample_umap
from .annotation import annotate_cells, filter_high_confidence, train_classifier, predict_all_cells
from .de_analysis import run_pseudobulk_de
from .gsea import run_gsea
from .io import load_adata

import scanpy as sc
import pandas as pd
from scipy.stats import mannwhitneyu
import numpy as np
import gc

class PostAnalysisPipeline:
    def __init__(self, h5ad_path, model_path):
        self.h5ad_path = h5ad_path
        self.model_path = model_path
        self.output_dir = ensure_dir(DEFAULT_ANALYSIS_ROOT / Path(model_path).parents[0].name)

    def paths(self, name):
        return self.output_dir / FILES[name]

    def run_latent(self, batch_size=1024):
        print("run latent...")
        return compute_scvi_latent(
            self.h5ad_path,
            self.model_path,
            self.paths('latent'),
            batch_size,
        )

    def run_clustering(self, resolution=1.0):
        print("run clustering...")
        return compute_leiden(
            self.paths('latent'),
            self.paths('leiden'),
            resolution,
        )

    def run_umap(
            self,
            output_path,
            adata_path=None,
            group_key="condition",
            groups=None,
            n_cells=100000,
            random_state=0,
            basis="X_scVI",
    ):
        """
        Génère un UMAP PNG.

        Parameters
        ----------
        output_path : str | Path
            Chemin du PNG de sortie.
        adata_path : str | Path | None
            h5ad source. None = self.paths("latent")
        group_key : str
            Colonne obs utilisée pour la couleur.
        groups : list | str | None
            Valeurs à afficher.
        n_cells : int
            Nombre max de cellules.
        random_state : int
            Seed.
        basis : str
            Embedding à utiliser.
        """
        print("run umap...")

        if adata_path is None:
            adata_path = self.paths("latent")

        return subsample_umap(
            latent_path=adata_path,
            output_png=output_path,
            group_key=group_key,
            groups=groups,
            n_cells=n_cells,
            random_state=random_state,
            basis=basis,
        )

    def run_annotation(self, marker_dict):
        print("run annotation...")
        print("annotate cells...")
        annotate_cells(
            self.paths('leiden'),
            self.paths('annotated'),
            marker_dict,
        )
        gc.collect()
        print("filter high confidence...")
        filter_high_confidence(
            self.paths('annotated'),
            self.paths('high_conf'),
        )
        gc.collect()
        print("train classifier...")
        train_classifier(
            self.paths('high_conf'),
            self.output_dir / 'celltype_classifier.pkl',
        )
        gc.collect()
        predict_all_cells(
            self.paths("latent"),
            self.output_dir / 'celltype_classifier.pkl',
            self.output_dir / "adata_final_annotated.h5ad",
        )


    def run_de(self, split_by=('cell_type', 'region')):
        print("run de...")
        return run_pseudobulk_de(
            self.paths('latent'),
            self.output_dir / 'celltype_classifier.pkl',
            self.output_dir / 'DE_results',
            split_by=split_by,
        )

    def run_gsea(self):
        print("run gsea...")
        return run_gsea(
            self.output_dir / 'DE_results',
            self.output_dir / 'GSEA_results',
        )

    def run_subcluster_analysis(
            self,
            target_cell_type="Inhibitory",
            resolution=0.5,
            min_cells_cluster=20,
            split_by=("subcluster",)
    ):
        """
        Subclustering d'un cell type puis DE pseudobulk par
        (subcluster x region), GSEA + figures UMAP.

        Tous les résultats sont sauvegardés dans :
        output/subclusters/<target_cell_type>/
        """

        print(f"\n=== Subcluster analysis: {target_cell_type} ===")

        # --------------------------------------------------
        # Paths
        # --------------------------------------------------
        root = self.output_dir / "subclusters"
        cell_dir = root / target_cell_type

        # reset uniquement ce dossier spécifique
        reset_folder(cell_dir)

        de_dir = cell_dir / "DE"
        gsea_dir = cell_dir / "GSEA"
        fig_dir = cell_dir / "figures"
        temp_h5ad = cell_dir / f"{target_cell_type}_subclusters.h5ad"

        de_dir.mkdir(parents=True, exist_ok=True)
        gsea_dir.mkdir(parents=True, exist_ok=True)
        fig_dir.mkdir(parents=True, exist_ok=True)

        # --------------------------------------------------
        # Load full data
        # --------------------------------------------------
        adata = load_adata(self.paths("annotated"))

        if "cell_type" not in adata.obs.columns:
            raise ValueError("cell_type column missing")

        # --------------------------------------------------
        # Filter target cell type
        # --------------------------------------------------
        ad = adata[adata.obs["cell_type"] == target_cell_type].copy()

        if ad.n_obs == 0:
            raise ValueError(f"No cells found for {target_cell_type}")

        print(f"{target_cell_type}: {ad.n_obs} cells")

        # --------------------------------------------------
        # Local clustering
        # --------------------------------------------------
        sc.pp.neighbors(ad, use_rep="X_scVI")
        sc.tl.leiden(ad, resolution=resolution, key_added="subcluster")

        ad.obs["subcluster"] = ad.obs["subcluster"].astype(str)

        # remove tiny clusters
        vc = ad.obs["subcluster"].value_counts()
        keep = vc[vc >= min_cells_cluster].index
        ad = ad[ad.obs["subcluster"].isin(keep)].copy()

        print("Kept clusters:")
        print(ad.obs["subcluster"].value_counts())

        # --------------------------------------------------
        # Cluster proportions ALS vs CTRL per patient
        # --------------------------------------------------
        print("\nComputing cluster proportions per patient...")

        required = ["patient_id", "condition", "subcluster"]
        for col in required:
            if col not in ad.obs.columns:
                raise ValueError(f"{col} missing in obs")

        tmp = ad.obs[required].copy()

        totals = tmp.groupby("patient_id").size().rename("total_cells")

        counts = (
            tmp.groupby(["patient_id", "condition", "subcluster"])
            .size()
            .rename("n_cells")
            .reset_index()
        )

        counts = counts.merge(totals, on="patient_id")
        counts["proportion"] = counts["n_cells"] / counts["total_cells"]

        counts.to_csv(
            cell_dir / "cluster_proportions_per_patient.csv",
            index=False
        )

        rows = []
        for cl in counts["subcluster"].unique():
            sub = counts[counts["subcluster"] == cl]

            als = sub[sub["condition"] == "ALS"]["proportion"].values
            ctrl = sub[sub["condition"] == "CTRL"]["proportion"].values

            if len(als) >= 2 and len(ctrl) >= 2:
                stat, p = mannwhitneyu(als, ctrl, alternative="two-sided")
            else:
                p = None

            rows.append({
                "subcluster": cl,
                "ALS_mean": als.mean() if len(als) else 0,
                "CTRL_mean": ctrl.mean() if len(ctrl) else 0,
                "log2FC": np.log2((als.mean() + 1e-9) / (ctrl.mean() + 1e-9)),
                "pval": p
            })

        pd.DataFrame(rows).to_csv(
            cell_dir / "cluster_enrichment_ALS_vs_CTRL.csv",
            index=False
        )

        # --------------------------------------------------
        # Marker genes
        # --------------------------------------------------
        print("\nComputing marker genes...")

        ad_mark = ad.copy()
        sc.pp.normalize_total(ad_mark, target_sum=1e4)
        sc.pp.log1p(ad_mark)

        sc.tl.rank_genes_groups(
            ad_mark,
            groupby="subcluster",
            method="wilcoxon",
            pts=True
        )

        rg = ad_mark.uns["rank_genes_groups"]
        groups = rg["names"].dtype.names

        rows = []
        n_top = 20

        for g in groups:
            genes = rg["names"][g][:n_top]
            scores = rg["scores"][g][:n_top]
            pvals = rg["pvals_adj"][g][:n_top] if "pvals_adj" in rg else [None] * len(genes)
            lfc = rg["logfoldchanges"][g][:n_top] if "logfoldchanges" in rg else [None] * len(genes)

            for gene, score, padj, fc in zip(genes, scores, pvals, lfc):
                rows.append({
                    "subcluster": g,
                    "gene": gene,
                    "score": score,
                    "logFC": fc,
                    "pval_adj": padj
                })

        pd.DataFrame(rows).to_csv(
            cell_dir / "top_markers.csv",
            index=False
        )

        # --------------------------------------------------
        # Save temp adata
        # --------------------------------------------------
        ad.write(temp_h5ad)

        # --------------------------------------------------
        # UMAP figures
        # --------------------------------------------------
        print("\nSaving UMAP figures...")

        for key in ["condition", "subcluster", "region"]:
            if key in ad.obs.columns:
                self.run_umap(
                    output_path=fig_dir / f"umap_{key}.png",
                    adata_path=temp_h5ad,
                    group_key=key,
                    n_cells=100000,
                )

        # --------------------------------------------------
        # Run DE
        # --------------------------------------------------
        print("\nRunning pseudobulk DE...")

        run_pseudobulk_de(
            latent_path=temp_h5ad,
            classifier_path=self.output_dir / "celltype_classifier.pkl",
            output_dir=de_dir,
            split_by=split_by,
            min_cells=3,
            min_total_cells=10,
            min_patients=2,
            reclassify=False
        )

        # --------------------------------------------------
        # Run GSEA
        # --------------------------------------------------
        print("\nRunning GSEA...")

        run_gsea(
            de_dir=de_dir,
            output_dir=gsea_dir,
        )

        print("\nDone.")
        print(cell_dir)