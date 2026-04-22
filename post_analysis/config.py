from pathlib import Path
DEFAULT_ANALYSIS_ROOT = Path('./analysis')
FILES = {
 'latent':'scvi_latent.h5ad',
 'leiden':'adata_with_leiden.h5ad',
 'annotated':'adata_annotated.h5ad',
 'high_conf':'adata_high_confidence.h5ad'
}