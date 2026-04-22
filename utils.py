import os
import shutil
def get_next_run_dir(base_dir="./runs", prefix="run"):
    os.makedirs(base_dir, exist_ok=True)

    existing = []

    for name in os.listdir(base_dir):
        if name.startswith(prefix + "_"):
            try:
                idx = int(name.split("_")[1])
                existing.append(idx)
            except:
                pass

    i = 0
    while i in existing:
        i += 1

    run_dir = os.path.join(base_dir, f"{prefix}_{i}")
    os.makedirs(run_dir, exist_ok=False)  # crée le dossier

    return run_dir

def reset_folder(path):
    if os.path.exists(path):
        shutil.rmtree(path)  # supprime tout le dossier
    os.makedirs(path)