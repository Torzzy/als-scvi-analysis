from pathlib import Path
import shutil

def ensure_dir(path):
    p=Path(path); p.mkdir(parents=True, exist_ok=True); return p

def reset_folder(path):
    p=Path(path)
    if p.exists(): shutil.rmtree(p)
    p.mkdir(parents=True, exist_ok=True)
    return p