from typing import Dict, List, Set
import pandas as pd
from collections import defaultdict
import numpy as np

ALS_REFERENCE = {

    # METABOLISM
    "Oxidative Phosphorylation": {
        "ATP5F1A","ATP5F1B","ATP5MC3","CYCS","UQCRB","UQCRH",
        "NDUFA4","NDUFA5","NDUFS1","COX4I1","COX7A2L",
        "VDAC1","SLC25A3","SLC25A4","MDH1","LDHB","ISCU"
    },

    "TCA Cycle": {
        "MDH1","IDH3A","SUCLA2","PDHA1","PDHB","CS","ACO2","DLST"
    },

    "Respiratory Electron Transport": {
        "ATP5F1A","ATP5F1B","CYCS","UQCRB","UQCRH","COX4I1",
        "NDUFA4","NDUFS1","COX7A2L","ATP5MC3"
    },

    # STRESS / PROTEOSTASIS
    "Heat Shock Response": {
        "HSPA1A","HSPA8","HSPA5","HSP90AA1","HSP90AB1",
        "HSPH1","HSPB1","DNAJB1","PTGES3","HSPA4L"
    },

    "Unfolded Protein Response": {
        "HSPA5","ERN1","EIF2AK3","ATF6","DDIT3",
        "XBP1","DNAJB9","HERPUD1","CANX"
    },

    # RNA / GENETIC REGULATION (ALS core)
    "RNA processing": {
        "HNRNPA1","HNRNPA2B1","HNRNPU","SRSF1","SRSF3",
        "DDX3X","DDX5","TARDBP","FUS","RBM14","TIA1"
    },

    "mRNA splicing": {
        "SRSF1","SRSF3","SF3B1","U2AF2","HNRNPA1",
        "HNRNPU","PRPF8","DDX5","RBM25"
    },

    "MYC targets": {
        "EIF4G2","EIF1AX","RAN","XPO1","PGK1","GAPDH",
        "RPLP0","RPS3","HNRNPA3","NPM1","MYC"
    },

    # NEURONAL FUNCTION
    "Glutamate signaling": {
        "SLC17A7","SLC17A6","GRIN1","GRIA1","GRIN2B",
        "SLC1A2","SLC1A3","GLS","GLUL"
    },

    "Synaptic transmission": {
        "SNAP25","SYT1","CAMK2A","GABRA1","GAD1","GAD2",
        "DLG4","SLC17A7","SLC17A6","GRIN1"
    },

    # IMMUNE / INFLAMMATION
    "NF-kB signaling": {
        "NFKB1","RELA","NFKBIA","TNF","IL1B","IL6",
        "CXCL10","STAT1","TRAF6","IKBKB"
    },

    "Interferon response": {
        "STAT1","STAT2","IRF1","IRF7","IFIT1","IFIT3",
        "ISG15","OAS1","MX1","CXCL10"
    },

    "Antigen presentation": {
        "HLA-DRA","HLA-DRB1","HLA-DPA1","B2M",
        "TAP1","TAP2","CD74","CIITA"
    },

    # TRANSPORT / CYTOSKELETON
    "Axonal transport": {
        "KIF5A","KIF1A","DYNC1H1","DYNLL1","NEFL",
        "MAPT","RAB7A","RAB10","RAB14","TUBA1A"
    },

    "Vesicle trafficking": {
        "RAB5A","RAB7A","RAB10","RAB14","CLTC",
        "AP2M1","DNM1","SNAP25","STX1A"
    },


    # DEGRADATION SYSTEMS
    "Autophagy": {
        "MAP1LC3B","SQSTM1","ATG5","ATG7","BECN1",
        "ULK1","VPS34","GABARAP","OPTN"
    },

    "Proteasome": {
        "PSMA1","PSMA3","PSMB1","PSMB5","PSMC1",
        "PSMD1","UBC","UBB","VCP","CAND1"
    }
}




def jaccard(a: Set[str], b: Set[str]) -> float:
    if len(a) == 0 and len(b) == 0:
        return 0.0
    return len(a.intersection(b)) / len(a.union(b))


def score_celltype_vs_reference(
    lead_pathway_genes: Dict[str, List[str]],
    als_reference: Dict[str, List[str]]
) -> float:
    """
    Calcule un score d'enrichissement entre les pathways d’un celltype
    (issus des lead genes) et une référence ALS.

    :param lead_pathway_genes: dict {pathway -> liste de gènes} pour un celltype
    :param als_reference: dict {pathway -> liste de gènes} de référence ALS
    :return: score moyen des meilleurs Jaccard (float)
    """


    best_scores = []

    # conversion en sets (important)
    ref_sets = {
        k: set(v) for k, v in als_reference.items()
    }

    for lead_pw, lead_genes in lead_pathway_genes.items():

        lead_set = set(lead_genes)

        best_score = 0.0

        for ref_pw, ref_genes in ref_sets.items():

            score = jaccard(lead_set, ref_genes)

            if score > best_score:
                best_score = score

        best_scores.append(best_score)

    if len(best_scores) == 0:
        return 0.0

    return sum(best_scores) / len(best_scores)

def compute_celltype_scores_from_csv(
    csv_path: str,
    als_reference: dict,
    scoring_function,
    contrast_filter: str = "global_condition"
):
    """
    Calcule un score d’enrichissement ALS pour chaque celltype à partir
    d’un fichier CSV contenant des résultats de pathways (ex: GSEA).

    Étapes :
    - Charge le CSV et filtre sur un contraste donné
    - Regroupe les gènes "lead" par celltype et par pathway
    - Pour chaque celltype :
        - applique une fonction de scoring (ex: Jaccard vs référence ALS)
        - agrège tous les gènes (tous pathways confondus)

    :param csv_path: chemin vers le fichier CSV contenant les pathways enrichis
    :param als_reference: dict {pathway -> liste de gènes} de référence ALS
    :param scoring_function: fonction de scoring entre pathways et référence
    :param contrast_filter: valeur de la colonne "contrast" à filtrer
    :return: DataFrame avec colonnes [celltype, score, lead_genes]
    """

    df = pd.read_csv(csv_path)
    df = df[df["contrast"] == contrast_filter]

    celltype_to_pathways = defaultdict(dict)

    for _, row in df.iterrows():

        celltype = row["celltype"]
        pathway = row["pathway"]

        genes = str(row["lead_genes"]).split(";")
        genes = [g.strip() for g in genes if g.strip()]

        celltype_to_pathways[celltype][pathway] = genes

    results = []

    for celltype, lead_dict in celltype_to_pathways.items():

        score = scoring_function(
            lead_dict,
            als_reference
        )

        # flatten all lead genes for this celltype
        all_genes = []
        for gset in lead_dict.values():
            all_genes.extend(gset)

        results.append({
            "celltype": celltype,
            "score": score,
            "lead_genes": list(set(all_genes))
        })

    return pd.DataFrame(results)

def permutation_test_celltypes(
    celltype_df,
    als_reference,
    scoring_function,
    strata,
    gene_df,
    n_iter=1000,
    random_state=0
):
    """
    Effectue un test de permutation (Monte Carlo) pour chaque celltype
    afin d’estimer une distribution nulle de scores.

    Étapes :
    - Pour chaque celltype :
        - récupère les lead genes observés
        - génère des ensembles de gènes aléatoires appariés (stratifiés)
          selon les propriétés des lead genes (ex: expression moyenne)
        - calcule un score pour chaque tirage aléatoire
    - Estime en ligne (streaming) la moyenne et la variance des scores nuls
      sans stocker toute la distribution

    :param celltype_df: DataFrame avec colonnes ["celltype", "Lead_genes"]
    :param als_reference: dict {pathway -> liste de gènes} de référence ALS
    :param scoring_function: fonction de scoring entre gènes et référence
    :param strata: dict {bin -> liste de gènes} pour le sampling stratifié
    :param gene_df: DataFrame contenant au moins ["gene", "bin"]
    :param n_iter: nombre d’itérations de permutation
    :param random_state: seed pour reproductibilité
    :return: dict {celltype -> {"mean", "std", "scores", "k"}}
    """


    rng = np.random.default_rng(random_state)

    results = {}

    for _, row in celltype_df.iterrows():
        print(f"Celltype : {row['celltype']}")
        celltype = row["celltype"]
        lead_genes = row["Lead_genes"]

        k = len(lead_genes)

        # streaming stats
        mean = 0.0
        m2 = 0.0
        count = 0

        obs_scores = []

        # OBSERVE: optional external pass
        for _ in range(n_iter):
            sampled_genes = matched_stratified_sample(
                strata,
                lead_genes,
                gene_df,
                rng
            )

            fake_lead = {"random": list(sampled_genes)}
            score = scoring_function(fake_lead, als_reference)

            obs_scores.append(score)

            count += 1
            delta = score - mean
            mean += delta / count
            m2 += delta * (score - mean)

        var = m2 / (count - 1) if count > 1 else 0
        std = np.sqrt(var)

        results[celltype] = {
            "mean": mean,
            "std": std,
            "scores": obs_scores,  # optional only if small
            "k": k
        }

    return results

def build_gene_strata_with_features(adata, n_bins=7):
    """
    Construit des strates de gènes basées sur leur niveau d’expression moyen
    afin de permettre un échantillonnage stratifié.

    Étapes :
    - Calcule l’expression moyenne de chaque gène sur toutes les cellules
    - Regroupe les gènes en quantiles (bins) d’expression
    - Assigne chaque gène à une strate correspondant à son bin

    :param adata: objet AnnData contenant les données d’expression (cells x genes)
    :param n_bins: nombre de bins (quantiles) pour discrétiser l’expression
    :return:
        - strata: dict {bin -> liste de gènes}
        - df: DataFrame avec colonnes ["gene", "mean_expr", "bin"]
    """
    gene_names = adata.var_names.to_numpy()

    mean_expr = np.asarray(adata.X.mean(axis=0)).ravel()

    df = pd.DataFrame({
        "gene": gene_names,
        "mean_expr": mean_expr
    })

    df["bin"] = pd.qcut(df["mean_expr"], q=n_bins, duplicates="drop")

    strata = defaultdict(list)

    for _, row in df.iterrows():
        strata[row["bin"]].append(row["gene"])

    return strata, df


def get_gene_signature(adata, genes):
    """
    Calcule des statistiques descriptives pour un ensemble de gènes
    à partir des données d’expression.

    Étapes :
    - Récupère les indices des gènes dans l’objet AnnData
    - Calcule l’expression moyenne de chaque gène sur toutes les cellules
    - Extrait les valeurs correspondant aux gènes d’intérêt
    - Résume ces valeurs via une moyenne et des quantiles

    :param adata: objet AnnData contenant les données d’expression (cells x genes)
    :param genes: liste de gènes pour lesquels calculer la signature
    :return: dict avec :
        - "mean_expr": moyenne des expressions des gènes
        - "expr_quantiles": quantiles (25%, 50%, 75%) des expressions
    """
    gene_idx = adata.var_names.get_indexer(genes)

    expr = np.asarray(adata.X.mean(axis=0)).ravel()

    return {
        "mean_expr": expr[gene_idx].mean(),
        "expr_quantiles": np.percentile(expr[gene_idx], [25, 50, 75])
    }


def matched_stratified_sample(strata, lead_genes, gene_df, rng):
    """
    Génère un échantillon de gènes aléatoires apparié aux lead genes
    en respectant leur distribution selon les strates (bins).

    Étapes :
    - Associe chaque lead gene à son bin d’expression
    - Pour chaque bin :
        - tire aléatoirement un gène dans la même strate
    - Retourne l’ensemble des gènes tirés (unicité assurée)

    :param strata: dict {bin -> liste de gènes}
    :param lead_genes: liste de gènes observés à apparier
    :param gene_df: DataFrame contenant au moins ["gene", "bin"]
    :param rng: générateur aléatoire numpy
    :return: liste de gènes échantillonnés
    """
    lead_bins = gene_df.set_index("gene").loc[lead_genes, "bin"]

    sampled = []

    for b in lead_bins:
        b_genes = strata[b]
        sampled.append(rng.choice(b_genes))

    return list(set(sampled))



