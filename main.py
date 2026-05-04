
import train_scvi
from post_analysis.pipeline import PostAnalysisPipeline
from preprocessing.pipeline_preprocessing import PipelinePreprocessing


if __name__ == '__main__':
    aux = 2
    if aux == 0:
        # fais le preprocessing, adata pret à l'entrainement du scvi.
        pp = PipelinePreprocessing(root_path='./root_data')
        pp.run()
    elif aux == 1:
        # entrainement du modèle
        train_scvi.train_scvi('./root_data/scvi_ready_data/final_scvi.h5ad',)
    elif aux == 2:
        # analyse
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

        pipe = PostAnalysisPipeline('./root_data/scvi_ready_data/final_scvi.h5ad', "./runs/run_3/model/")
        #pipe.run_latent()
        #pipe.run_clustering(resolution=1.0)
        #pipe.run_cell_annotation(MARKERS)
        #pipe.run_metrics(MARKERS)
        #pipe.run_de()
        pipe.run_gsea(output_dir='GSEA_relaxed_results', min_genes=5, fdr_threshold=.1, nes_threshold=1.2)