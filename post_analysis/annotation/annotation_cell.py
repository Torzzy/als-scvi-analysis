import numpy as np
import scanpy as sc
from pathlib import Path

from post_analysis.io import save_adata
from post_analysis.config import FILES

def annotate_cells_cluster_aware(
    adata,
    marker_dict,
    output_dir,
    cluster_key="leiden",
    latent_key="X_scVI",
    neighbors_key="connectivities",
    n_neighbors=15,
    max_iter=10,
    tol=1e-3,
    weights=(0.5, 0.3, 0.2),  # neigh, marker, centroid
):
    """
    Cluster-aware RAM-safe iterative annotation.

    Auto-builds neighbors if missing.

    Outputs:
        adata.obs["celltype_pred"]
        adata.obs["celltype_confidence"]
    """

    w_neigh, w_mark, w_cent = weights

    n_cells = adata.n_obs
    cell_types = list(marker_dict.keys())
    n_types = len(cell_types)

    # ======================================================
    # 0. BUILD NEIGHBORS IF MISSING
    # ======================================================
    if neighbors_key not in adata.obsp:
        print("Building neighbors graph...")
        sc.pp.neighbors(
            adata,
            use_rep=latent_key,
            n_neighbors=n_neighbors,
            key_added="neighbors"
        )

    conn = adata.obsp["connectivities"]

    # ======================================================
    # 1. INPUTS
    # ======================================================
    latent = adata.obsm[latent_key]
    clusters = adata.obs[cluster_key].values

    # ======================================================
    # 3. MARKER SCORES (RAM SAFE)
    # ======================================================
    print("Computing marker scores...")

    gene_to_idx = {g: i for i, g in enumerate(adata.var_names)}
    marker_scores = np.zeros((n_cells, n_types), dtype=np.float32)

    X = adata.X

    for j, ct in enumerate(cell_types):

        idx = [gene_to_idx[g] for g in marker_dict[ct] if g in gene_to_idx]

        if len(idx) == 0:
            continue

        X_sub = X[:, idx]

        if hasattr(X_sub, "toarray"):
            X_sub = X_sub.toarray()

        marker_scores[:, j] = np.asarray(X_sub.mean(axis=1)).ravel()

    marker_scores /= (marker_scores.max(axis=1, keepdims=True) + 1e-9)

    # ======================================================
    # 4. INIT SCORES
    # ======================================================
    scores = marker_scores.copy()

    unique_clusters = np.unique(clusters)

    # ======================================================
    # 5. ITERATIVE REFINEMENT
    # ======================================================
    for it in range(max_iter):

        print(f"Iteration {it+1}")

        # -------------------------
        # neighbor propagation
        # -------------------------
        neigh_scores = conn.dot(scores)

        row_sum = np.asarray(conn.sum(axis=1)).ravel()[:, None] + 1e-9
        neigh_scores /= row_sum

        # -------------------------
        # centroid per (cluster, type)
        # -------------------------
        best_idx = np.argmax(scores, axis=1)

        centroids = {}

        for cl in unique_clusters:
            idx_cl = np.where(clusters == cl)[0]

            for j in range(n_types):

                idx = idx_cl[best_idx[idx_cl] == j]

                if len(idx) == 0:
                    continue

                centroids[(cl, j)] = latent[idx].mean(axis=0)

        centroid_scores = np.zeros((n_cells, n_types), dtype=np.float32)

        for i in range(n_cells):
            cl = clusters[i]

            for j in range(n_types):

                key = (cl, j)

                if key not in centroids:
                    continue

                dist = np.linalg.norm(latent[i] - centroids[key])
                centroid_scores[i, j] = 1 / (1 + dist)

        centroid_scores /= (centroid_scores.max(axis=1, keepdims=True) + 1e-9)

        # -------------------------
        # combine
        # -------------------------
        new_scores = (
            w_neigh * neigh_scores +
            w_mark * marker_scores +
            w_cent * centroid_scores
        )

        new_scores /= (new_scores.max(axis=1, keepdims=True) + 1e-9)

        # -------------------------
        # convergence
        # -------------------------
        delta = np.abs(new_scores - scores).mean()
        print(f"delta = {delta:.6f}")

        scores = new_scores

        if delta < tol:
            print("Converged")
            break

    # ======================================================
    # 6. FINAL ASSIGNMENT
    # ======================================================
    best_idx = np.argmax(scores, axis=1)
    best_scores = scores[np.arange(n_cells), best_idx]

    sorted_scores = np.sort(scores, axis=1)
    second_scores = sorted_scores[:, -2]

    confidence = 1 / (1 + np.exp(-(best_scores - second_scores)))

    pred_labels = np.array(cell_types)[best_idx]

    # ======================================================
    # 7. STORE
    # ======================================================
    adata.obs["celltype_pred"] = pred_labels
    adata.obs["celltype_confidence"] = confidence

    # ======================================================
    # 8. SAVE
    # ======================================================
    save_path = Path(output_dir) / FILES["annotated"]
    save_adata(adata,save_path)

    return adata