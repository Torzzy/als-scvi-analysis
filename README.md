# Projet : Etude multi-cohortes de données snRNA-seq  via scVI pour une identification locale des mécanismes biologiques impliqués dans la SLA.
<img src="sources/pipe_pres.png" width="800">
<p align="center">
<b></b> Pipeline d’analyse snRNA-seq multi-cohortes basé sur scVI, edgeR et GSEA.
</p>

## English summary

This project analyzes single-nucleus RNA-seq datasets in Amyotrophic Lateral Sclerosis (ALS) using scVI for multi-cohort integration, cell-type annotation, differential expression, and pathway enrichment analysis.
The goal is to identify robust and tissue-specific transcriptional alterations across independent datasets.
## Résumé
Ce projet analyse des données de snRNA-seq dans la sclérose latérale amyotrophique (SLA) à l’aide de scVI pour l’intégration multi-cohortes, l’annotation des types cellulaires, l’analyse d’expression différentielle et l’enrichissement en voies biologiques.
L’objectif est d’identifier des altérations transcriptionnelles robustes et spécifiques aux tissus à travers plusieurs jeux de données indépendants.
## Contributions principales

- Intégration multi-cohortes de données snRNA-seq via scVI avec correction des batch effects
- Construction d’un espace latent biologique pour l’analyse cellulaire non supervisée
- Annotation hiérarchique des types cellulaires (markers + structure du latent space)
- Analyse différentielle en pseudobulk avec edgeR et modèle condition × région
- Enrichment de pathways (GSEA preranked) à résolution cellulaire
- Validation des signatures biologiques par comparaison à la littérature via tests permutationnels

## Objectif biologique

L’objectif de ce projet est de caractériser les altérations transcriptionnelles associées à la SLA dans différents types cellulaires et régions cérébrales (cortex moteur, cortex frontal, moelle épinière), à partir de données snRNA-seq multi-cohortes.

## Jeux de données

Trois jeux de données GEO ont été utilisés :

- GSE219280 : cortex frontal (FCX) et moteur (MCX), 12 individus
- GSE287257 : cortex préfrontal (FX), 12 individus
- GSE290359 : moelle épinière (SC), 48 individus

## Pipeline global

1. Prétraitement des données (QC, sélection des gènes, HVG)
2. Intégration multi-cohortes avec scVI
3. Analyse de l’espace latent et évaluation des batch effects
4. Clustering non supervisé (Leiden)
5. Annotation hiérarchique des cellules
6. Agrégation en pseudobulk par patient et type cellulaire
7. Analyse différentielle (edgeR)
8. Enrichment de pathways (GSEA preranked)
9. Validation des résultats par comparaison à la littérature

## Structure du projet

- `docs/methods.md` → description complète et détaillée de la méthode (niveau article scientifique)
- `docs/results.md` → résultats biologiques, métriques et figures
- `sources/` → figures, dotplots, tableaux de pathways

## Comment explorer ce projet

Pour une lecture rapide des principaux résultats :
- consulter `docs/results.md`

Pour une lecture complète :
- lire `docs/methods.md`

Sections importantes dans `docs/methods.md` :
- [données](docs/methods.md#donnees)
- [pré-traitement](docs/methods.md#pre-traitement-pour-scvi)
- [présentation scVI](docs/methods.md#scvi) 
- [annotation cluster/cellule](docs/methods.md#annotation-cellulaire) 
- [pseudobulk DE](docs/methods.md#pseudobulk-de) 
- [GSEA](docs/methods.md#gsea) 
- [Validationsur la littérature existante](docs/methods.md#validation)

## Résultats principaux

- Enrichissement significatif de pathways associés à la SLA dans les neurones et astrocytes
- Cohérence globale avec les signatures décrites dans la littérature

## Limites

- Annotation cellulaire dépendante de marqueurs et d’approches heuristiques (absence d’atlas de référence)
- Risque de mélange entre effets biologiques régionaux et batch effects résiduels malgré scVI
- Interprétation des pathways assistée par LLM, nécessitant prudence biologique

## Get started

Voici un exemple de comment utiliser le code :
```python
from preprocessing.pipeline_preprocessing import PipelinePreprocessing
pp = PipelinePreprocessing(root_path="./root_data")
pp.run()

import train_scvi
train_scvi.train_scvi("./root_data/scvi_ready_data/final_scvi.h5ad")

from post_analysis.pipeline import PostAnalysisPipeline
pipe = PostAnalysisPipeline("./root_data/scvi_ready_data/final_scvi.h5ad", "./runs/run_n/model/")
pipe.run_latent()
pipe.run_clustering(resolution=1.0)
MARKERS = {
    "Neuron": ["RBFOX3", "SYT1", "SNAP25", "MAP2"],
    "Excitatory": ["SLC17A7", "CAMK2A"],
    "Inhibitory": ["GAD1", "GAD2"],
    "Astrocyte": ["GFAP", "AQP4", "SLC1A2"],
    "Microglia": ["C1QA", "C1QB", "P2RY12", "TMEM119"],
    "Oligodendrocyte": ["MBP", "MOG", "PLP1"],
    "OPC": ["PDGFRA", "CSPG4"],
    "Endothelial": ["CLDN5", "FLT1", "PECAM1"],
    "Pericyte": ["PDGFRB", "RGS5"]
}

pipe.run_cluster_annotation(MARKERS)
pipe.run_metrics()
pipe.run_de()
pipe.run_pathway_validation()
```

Pour voir comment configurer le projet et lancer le code : [Get started](docs/methods.md#get-started) 
## Documentation complète

-  Méthodes détaillées : [methods.md](docs/methods.md)
-  Résultats complets : [results.md](docs/results.md)
-  Explication de edgeR : [edgeR.md](docs/edgeR.md)
## Ressources matérielles

Le pipeline peut être memory-intensive, notamment lors de l'intégration scVI et des étapes de pseudobulk.  
Il est recommandé d’exécuter les étapes de manière séquentielle.

Le pipeline a été développé et testé sur :

* **OS** : Ubuntu 24.04
* **RAM** : 16 Go
* **GPU** : NVIDIA GeForce RTX 5060 (8 Go VRAM)
* **CPU** : AMD Ryzen 5 8400F

Pour toute question ou remarque :

- GitHub : @Torzzy
- Email : tomdauve@gmail.com