

from .config import DEFAULT_ANALYSIS_ROOT, FILES
from .utils import ensure_dir, reset_folder
from .latent import compute_scvi_latent
from .clustering import compute_leiden, subsample_umap
from .de_analysis import run_pseudobulk_de, summarize_edger_results
from .gsea import run_gsea
from .io import load_adata
from .annotation.annotation_cluster import build_pseudobulk_matrix, run_cluster_vs_rest_de_fast, annotate_clusters_from_de
from .annotation.annotation_cell import annotate_cells_cluster_aware
from .annotation.metrics import (compute_knn_purity, compute_silhouette_per_cluster, compute_cluster_entropy,
                                 compute_patient_mixing, compute_dotplot_scores, plot_dotplot_from_scores, compute_auc_separation)
from .pathway_validation import score_celltype_vs_reference, ALS_REFERENCE, permutation_test_celltypes, build_gene_strata_with_features

from pathlib import Path
import scanpy as sc
import pandas as pd
from scipy.stats import mannwhitneyu
import numpy as np
import os
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
            adata_path="latent",
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
        else:
            adata_path = self.paths(adata_path)

        subsample_umap(
            latent_path=adata_path,
            output_png=output_path,
            group_key=group_key,
            groups=groups,
            n_cells=n_cells,
            random_state=random_state,
            basis=basis,
        )

    def run_cluster_annotation(self, marker_dict, cluster_key="leiden"):
        """
        Full cluster annotation pipeline:
        - cluster vs rest DE (RAM-safe)
        - marker scoring
        - cluster labeling
        - save all results in output_dir/annotations
        """

        print("\n=== Running cluster annotation pipeline ===")

        annot_dir = self.output_dir / "annotations"
        de_dir = annot_dir / "de_results"

        annot_dir.mkdir(parents=True, exist_ok=True)
        de_dir.mkdir(parents=True, exist_ok=True)

        # --------------------------------------------------
        # 0. Load clustered data (IMPORTANT FIX)
        # --------------------------------------------------
        adata = load_adata(self.paths("leiden"))

        if cluster_key not in adata.obs.columns:
            raise ValueError(f"Missing cluster_key: {cluster_key}")

        # --------------------------------------------------
        # 1. RAM-safe pseudobulk precomputation (CRITICAL FIX)
        # --------------------------------------------------
        print("Building pseudobulk matrix...")



        X_pb, meta_pb, genes = build_pseudobulk_matrix(
            adata,
            groupby=(cluster_key, "region"),
            patient_key="patient_id"
        )

        # attach cluster column explicitly for DE step
        meta_pb["cluster"] = meta_pb[cluster_key]

        # --------------------------------------------------
        # 2. Cluster vs rest DE (NO adata slicing anymore)
        # --------------------------------------------------
        print("Running cluster vs rest DE...")


        run_cluster_vs_rest_de_fast(
            X_pb=X_pb,
           meta=meta_pb,
            genes=genes,
            output_dir=str(de_dir),
            cluster_key="cluster"
        )

        # --------------------------------------------------
        # 3. Marker-based annotation
        # --------------------------------------------------
        print("Computing marker-based annotation...")

        annotations = annotate_clusters_from_de(
            de_folder=str(de_dir),
            marker_dict=marker_dict,
        )

        # --------------------------------------------------
        # 4. Save outputs
        # --------------------------------------------------
        annotations.to_csv(
            annot_dir / "final_cluster_labels.csv",
            index=False
        )

        cluster_map = dict(zip(
            annotations["cluster"],
            annotations["assigned_type"]
        ))

        pd.Series(cluster_map).to_json(
            annot_dir / "cluster_to_celltype.json"
        )

        print(f"Done → {annot_dir}")

        return annotations

    def run_cell_annotation(self, marker_dict):
        adata = load_adata(self.paths("leiden"))
        annotate_cells_cluster_aware(adata,
            marker_dict,
            cluster_key="leiden",
            latent_key="X_scVI",
            output_dir=self.output_dir,
            neighbors_key="connectivities",
            n_neighbors=15,
            max_iter=10,
            tol=1e-3,
            weights=(0.5, 0.3, 0.2),  # neigh, marker, centroid
            )

    def run_metrics(self, marker_dict):
        path_metric_annotations = Path(self.output_dir) / "annotations" / "metrics"
        os.makedirs(path_metric_annotations, exist_ok=True)

        adata = load_adata(self.paths("annotated"))
        mapping = pd.read_csv(self.output_dir / "annotations"/ "final_cluster_labels.csv")

        cluster_to_type = dict(zip(
            mapping["cluster"].astype(str),
            mapping["assigned_type"]
        ))

        adata.obs["assigned_type_from_cluster"] = adata.obs["leiden"].astype(str).map(cluster_to_type)

        print("silhouette...")
        sil = compute_silhouette_per_cluster(adata)
        sil.to_csv(f"{path_metric_annotations}/silhouette.csv")

        print("AUC separation...")
        auc_sep = compute_auc_separation(adata, marker_dict)
        auc_sep.to_csv(f"{path_metric_annotations}/auc_separation.csv", index=False)

        print("knn purity...")
        knn = compute_knn_purity(adata)
        knn.to_csv(f"{path_metric_annotations}/knn_purity.csv")

        print("patient mixing...")
        mix = compute_patient_mixing(adata)
        mix.to_csv(f"{path_metric_annotations}/patient_mixing.csv")

        print("dotplot...")
        dot = compute_dotplot_scores(adata, marker_dict)
        dot.to_csv(f"{path_metric_annotations}/dotplot.csv", index=False)
        plot_dotplot_from_scores(dot, f"{path_metric_annotations}/dotplot.png")

        print("entropy...")
        ent = compute_cluster_entropy(adata)
        ent.to_csv(f"{path_metric_annotations}/entropy.csv", index=False)

        print("umap...")
        del adata
        gc.collect()
        self.run_umap(output_path=str(path_metric_annotations / "umap_by_celltype.png"),
                      adata_path='annotated',
                      group_key='celltype_pred')

        self.run_umap(output_path=str(path_metric_annotations / "umap_by_condition.png"),
                      adata_path='annotated',
                      group_key='condition')

        self.run_umap(output_path=str(path_metric_annotations / "umap_by_region.png"),
                      adata_path='annotated',
                      group_key='region')
        return {
            "knn": knn,
            "silhouette": sil,
            "mixing": mix,
            "dotplot": dot,
            "entropy": ent
        }


    def run_de(self):
        print("run de...")
        de_celltype_dir = Path(self.output_dir) / "DE_results"
        csv_path = os.path.join(self.output_dir,"DE_results",'merged_DE.csv')
        run_pseudobulk_de(
            self.paths('annotated'),
            de_celltype_dir,
        )
        summarize_edger_results(de_celltype_dir, csv_path)


    def run_gsea(self, output_dir='GSEA_results', min_genes=10, fdr_threshold=0.05, nes_threshold=1.5):
        print("run gsea...")
        return run_gsea(
            os.path.join(self.output_dir / 'DE_results', 'merged_DE.csv'),
            self.output_dir / output_dir,
            min_genes=min_genes,
            fdr_threshold=fdr_threshold,
            nes_threshold=nes_threshold,
        )

    def run_pathway_validation(self, contrast_filter='global_condition',
                               scoring_function=score_celltype_vs_reference,
                               n_iter=1000,
                               random_state=42
                               ):
        print("run pathway validation...")
        # PATHS
        csv_path = os.path.join(self.output_dir, "GSEA_results", "TOP_PATHWAYS.csv")
        als_reference = ALS_REFERENCE
        adata = load_adata(self.paths('annotated'))
        output_csv = Path(self.output_dir) / "pathway_validation.csv"

        strata, gene_df = build_gene_strata_with_features(adata, 7)

        # Load CSV input

        df = pd.read_csv(csv_path)
        df = df[df["contrast"] == contrast_filter]

        celltype_to_pathways = {}
        celltype_to_genes = {}


        for _, row in df.iterrows():

            ct = row["celltype"]
            pw = row["Term"]

            genes = str(row["Lead_genes"]).split(";")
            genes = [g.strip() for g in genes if g.strip()]

            if ct not in celltype_to_pathways:
                celltype_to_pathways[ct] = {}

            celltype_to_pathways[ct][pw] = genes

        # flatten genes per celltype
        for ct, pw_dict in celltype_to_pathways.items():
            all_genes = []
            for gset in pw_dict.values():
                all_genes.extend(gset)
            celltype_to_genes[ct] = list(set(all_genes))


        # OBSERVED SCORES
        observed_scores = {}

        for ct, lead_dict in celltype_to_pathways.items():
            score = scoring_function(lead_dict, als_reference)

            observed_scores[ct] = {
                "score": score,
                "k": len(celltype_to_genes[ct])
            }

        # NULL DISTRIBUTION

        null_results = permutation_test_celltypes(
            celltype_df=pd.DataFrame({
                "celltype": list(celltype_to_genes.keys()),
                "Lead_genes": list(celltype_to_genes.values())
            }),
            als_reference=als_reference,
            scoring_function=scoring_function,
            strata=strata,
            gene_df=gene_df,
            n_iter=n_iter,
            random_state=random_state
        )

        # BUILD FINAL TABLE
        rows = []

        for ct in celltype_to_genes.keys():
            print(null_results.keys())
            print(null_results[list(null_results.keys())[0]])
            obs = observed_scores[ct]["score"]
            null_scores = np.array(null_results[ct]["scores"])

            null_mean = null_scores.mean()
            null_std = null_scores.std(ddof=1)

            zscore = 0.0
            if null_std > 0:
                zscore = (obs - null_mean) / null_std

            p_emp = (np.sum(null_scores >= obs) + 1) / (len(null_scores) + 1)
            rows.append({
                "celltype": ct,
                "observed_score": obs,
                "null_mean": null_mean,
                "null_std": null_std,
                "zscore": zscore,
                "p_empirical": p_emp,
                "n_genes": observed_scores[ct]["k"]
            })

        result_df = pd.DataFrame(rows)

        # SAVE
        result_df.to_csv(output_csv, index=False)

        return result_df
