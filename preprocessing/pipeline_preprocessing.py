import os


import preprocessing.preprocessing_219280
import preprocessing.preprocessing_287257
import preprocessing.preprocessing_290359
import preprocessing.preprocessing_merge_for_scvi



class PipelinePreprocessing:
    """
    cette classe permet de preparer le .h5ad qui sera utilisé pour entrainer scvi.
    On suppose que le root_path est organisé de la manière suivante :
    root_path/data/219280
                  /287257
                  /290359
    avec chaque dossier qui correspond au dataset GEO associé. Plus précisément :
    - 219280 : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE219280, il faut télécharger les 34 cellBender_corrected_filtered.h5 des 34 samples du dataset et les mettre dans ce dossier. (pas besoin de prendre forcément les patients FTD)
    - 287257 : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287257, il faut télécharger les 12 filtered_feature_bc_matrix.h5 des 12 samples du datasetet les mettre dans ce dossier.
    - 290359 : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290359, il faut télécharger le fichier GSE290359_BA44_BA46_COUNTS.h5ad et le mettre dans le dossier.
    La classe va automatiquement générer des dossiers intermédiaires dans path_data. Le fichier .h5ad pret à l'emploi pour scvi se trouvera dans  root_path/scvi_ready_data/.
    """
    def __init__(self, root_path):
        self.root_path = root_path

    def prepare_data(self):
        """
        Permet d'uniformiser la structure des données pour pouvoir faire un traitement commun sur tous les fichiers .h5ad par la suite.
        :return: les données uniformisées pretes à etre utilisées par la fonction keep_common_genes dans root_path/preprocessed_data.
        """
        #paths
        data_path = os.path.join(self.root_path, "data")
        data_path_290359 = os.path.join(os.path.join(data_path, "290359"),'GSE290359_BA44_BA46_COUNTS.h5ad')
        data_path_287257 = os.path.join(data_path, "287257")
        data_path_219280 = os.path.join(data_path, "219280")

        # preparer chaque dataset independament
        preprocessing.preprocessing_290359.preprocess_all_patients(data_path_290359, self.root_path)
        preprocessing.preprocessing_287257.batch_preprocess_scvi_spinal(
            input_dir=data_path_287257,
            output_root=os.path.join(self.root_path,'preprocessed_data')
        )
        preprocessing.preprocessing_219280.batch_preprocess_scvi_auto(data_path_219280,os.path.join(self.root_path,'preprocessed_data'))

    def keep_common_genes(self):
        """
        sauvegarde les données en ne gardant que les genes communs à tous les fichiers.
        :return: les memes données filtrées dans root_path/common_genes_data
        """
        dataset_dirs = [
            os.path.join(self.root_path,"preprocessed_data/219280/"),
            os.path.join(self.root_path,"preprocessed_data/287257/"),
            os.path.join(self.root_path,"preprocessed_data/290359/")

        ]

        adata = preprocessing.preprocessing_merge_for_scvi.build_common_gene_datasets(
            dataset_dirs,
            output_root=os.path.join(self.root_path,"common_genes_data")
        )
    def perform_global_qc_hvg(self):
        """
        On calcule la moyenne et variance de chaque gene sur tous les datasets. on les tri par variances et on ne garde que les 3000 genes les plus variables.
        On ne garde que les cellules qui ont un minimum de genes exprimés (ici 200).
        On ne garde que les cellules qui ont moins de 20% de MT (pas stressées etc...)
        On vérifie qu'il n'y a pas de genes très rares présents seulement dans 1 cellule. (ici il faut qu'il y ait au moins 3 cellules qui expriment le gène pour qu'il soit gardé).

        :return: meme format que les données de départ stockées dans root_path/qc_hvg_data.
        """
        adata = preprocessing.preprocessing_merge_for_scvi.qc_hvg_global_from_root(
            os.path.join(self.root_path,"common_genes_data"),
            output_root=os.path.join(self.root_path, 'qc_hvg_data')
        )

    def merge_data(self):
        """
        On fusionne tous les .h5ad en un seul fichier directement utilisable par scvi.
        :return: le .h5ad final se trouve à root_path/scvi_ready_data/final_scvi.h5ad
        """
        os.makedirs(os.path.join(self.root_path,"scvi_ready_data"), exist_ok=True)
        adata = preprocessing.preprocessing_merge_for_scvi.merge_scvi_ready_datasets(
            os.path.join(self.root_path,"qc_hvg_data"),
            output_path=os.path.join(self.root_path,"scvi_ready_data/final_scvi.h5ad")
        )

    def run(self):
        self.prepare_data()
        self.keep_common_genes()
        self.perform_global_qc_hvg()
        self.merge_data()


