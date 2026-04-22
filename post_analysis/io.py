import scanpy as sc
from .utils import ensure_dir

def load_adata(path):
    return sc.read_h5ad(path)

def save_adata(adata, path, compression='gzip'):
    ensure_dir(__import__('pathlib').Path(path).parent)
    adata.write(path, compression=compression)