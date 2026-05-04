import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from sklearn.metrics import silhouette_samples
from sklearn.neighbors import NearestNeighbors

def compute_knn_purity(
    adata,
    label_col="celltype_pred",
    n_neighbors=15,
):
    X = adata.obsm["X_scVI"]

    nn = NearestNeighbors(n_neighbors=n_neighbors).fit(X)
    idx = nn.kneighbors(return_distance=False)

    labels = adata.obs[label_col].values

    purity = np.zeros(adata.n_obs, dtype=np.float32)

    for i in range(adata.n_obs):
        neigh_labels = labels[idx[i]]
        purity[i] = np.mean(neigh_labels == labels[i])

    adata.obs["knn_purity"] = purity

    return adata.obs.groupby("leiden")["knn_purity"].mean()

def compute_silhouette_per_cluster(
    adata,
    label_col="leiden",
    obsm_key="X_scVI",
    sample_max=5000,
    cluster_col="leiden"
):
    """
    Computes silhouette score correctly and safely.
    """

    X = adata.obsm[obsm_key]
    labels = adata.obs[label_col].astype(str).values
    clusters = adata.obs[cluster_col].astype(str).values

    n = adata.n_obs

    # --------------------------------------------------
    # subsampling (SAFE)
    # --------------------------------------------------
    if n > sample_max:
        idx = np.random.choice(n, sample_max, replace=False)
        X_sub = X[idx]
        labels_sub = labels[idx]
        clusters_sub = clusters[idx]
    else:
        X_sub = X
        labels_sub = labels
        clusters_sub = clusters

    # --------------------------------------------------
    # silhouette
    # --------------------------------------------------
    sil = silhouette_samples(X_sub, labels_sub)

    df = pd.DataFrame({
        "cluster": clusters_sub,
        "silhouette": sil
    })

    return df.groupby("cluster")["silhouette"].mean().reset_index()


import numpy as np
import pandas as pd
from scipy import sparse

def compute_auc_separation(
    adata,
    marker_dict,
    cluster_col="leiden",
    top_n=1000
):
    """
    True AUCell-like computation (rank-based).

    Steps:
    - rank genes per cell
    - compute enrichment of marker genes in top-ranked genes
    - aggregate per cluster
    """

    clusters = np.array(adata.obs[cluster_col].values)
    gene_index = {g: i for i, g in enumerate(adata.var_names)}
    X = adata.X

    results = []

    # ======================================================
    # iterate clusters (RAM-safe)
    # ======================================================
    for cl in np.unique(clusters):

        idx = np.where(clusters == cl)[0]
        if len(idx) == 0:
            continue

        Xc = X[idx]

        if sparse.issparse(Xc):
            Xc = Xc.toarray()

        n_cells = Xc.shape[0]
        n_genes = Xc.shape[1]

        # ======================================================
        # rank genes per cell
        # ======================================================
        ranked = np.argsort(-Xc, axis=1)  # descending

        # restrict to top genes for speed
        ranked = ranked[:, :top_n]

        type_scores = []

        # ======================================================
        # AUCell per cell type
        # ======================================================
        for cell_type, markers in marker_dict.items():

            marker_idx = [gene_index[g] for g in markers if g in gene_index]

            if len(marker_idx) == 0:
                type_scores.append(0.0)
                continue

            marker_set = set(marker_idx)

            # --------------------------------------------------
            # AUCell: fraction of marker genes in top ranks
            # --------------------------------------------------
            auc_values = []

            for i in range(n_cells):

                top_genes = ranked[i]

                overlap = len(set(top_genes) & marker_set)

                auc = overlap / len(marker_set)

                auc_values.append(auc)

            type_scores.append(np.mean(auc_values))

        type_scores = np.array(type_scores)

        if len(type_scores) < 2:
            continue

        best = np.max(type_scores)
        second = np.partition(type_scores, -2)[-2]

        results.append({
            "cluster": cl,
            "auc_separation": best - second,
            "auc_best": best,
            "auc_second": second
        })

    return pd.DataFrame(results)


def entropy(p):
    p = np.clip(p, 1e-9, 1)
    return -np.sum(p * np.log(p))


def compute_cluster_entropy(
    adata,
    conf_col="celltype_confidence",
    cluster_col="leiden"
):
    res = []

    for cl, sub in adata.obs.groupby(cluster_col):

        conf = sub[conf_col].values
        conf = conf / conf.sum()

        res.append({
            "cluster": cl,
            "entropy": entropy(conf)
        })

    return pd.DataFrame(res)

def compute_patient_mixing(
    adata,
    cluster_col="leiden",
    patient_col="patient_id",
    n_neighbors=15
):
    X = adata.obsm["X_scVI"]

    nn = NearestNeighbors(n_neighbors=n_neighbors).fit(X)
    idx = nn.kneighbors(return_distance=False)

    patients = adata.obs[patient_col].values
    clusters = adata.obs[cluster_col].values

    mixing = []

    for i in range(adata.n_obs):

        neigh_pat = patients[idx[i]]

        mix = len(set(neigh_pat)) / n_neighbors

        mixing.append(mix)

    adata.obs["patient_mixing"] = mixing

    return adata.obs.groupby(cluster_col)["patient_mixing"].mean()

def compute_dotplot_scores(
    adata,
    marker_dict,
    cluster_col="leiden"
):
    """
    RAM-safe dotplot:
    - no .toarray()
    - no per-cell loops
    """

    clusters = adata.obs[cluster_col].values
    X = adata.X  # keep sparse if possible

    gene_index = {g: i for i, g in enumerate(adata.var_names)}

    results = []

    for ct, markers in marker_dict.items():

        idx = [gene_index[g] for g in markers if g in gene_index]
        if len(idx) == 0:
            continue

        # extract only marker columns
        Xm = X[:, idx]

        if hasattr(Xm, "toarray"):  # sparse safe conversion ONLY on subset
            Xm = Xm.toarray()

        marker_score = Xm.mean(axis=1)

        df = pd.DataFrame({
            "cluster": clusters,
            "score": marker_score
        })

        agg = df.groupby("cluster")["score"].mean()

        for cl, sc in agg.items():
            results.append({
                "cluster": cl,
                "celltype": ct,
                "score": float(sc)
            })

    return pd.DataFrame(results)



def plot_dotplot_from_scores(
    dot_df,
    output_path,
    cluster_col="cluster",
    celltype_col="celltype",
    score_col="score"
):
    """
    Create dotplot PNG from compute_dotplot_scores output.

    Rows = clusters
    Columns = celltypes
    Values = mean marker score
    """

    # ----------------------------
    # pivot table
    # ----------------------------
    mat = dot_df.pivot(
        index=cluster_col,
        columns=celltype_col,
        values=score_col
    ).fillna(0)

    # ----------------------------
    # plot
    # ----------------------------
    plt.figure(figsize=(12, 6))

    im = plt.imshow(mat.values, aspect="auto")

    plt.xticks(
        range(len(mat.columns)),
        mat.columns,
        rotation=90
    )

    plt.yticks(
        range(len(mat.index)),
        mat.index
    )

    plt.colorbar(im, label="Marker score")

    plt.xlabel("Cell types")
    plt.ylabel("Clusters")
    plt.title("Dotplot (marker scores)")

    plt.tight_layout()

    plt.savefig(output_path, dpi=300)
    plt.close()

    return mat
