"""Mean number of enzyme and transporter genes per taxon in an input folder, over several random communities."""
import os, shutil, sys, tempfile

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
sys.path.insert(0, HERE)
from initialization import initialize_data  # noqa: E402
from sweep_minimal import set_runtime  # noqa: E402

src = os.path.abspath(sys.argv[1])
seeds = [int(s) for s in sys.argv[2:]] or [0, 1, 2]
enz, upt = [], []
for seed in seeds:
    tmp = tempfile.mkdtemp()
    w = os.path.join(tmp, "in")
    shutil.copytree(src, w)
    set_runtime(os.path.join(w, "runtime.txt"), dict(x=10, y=10, gridsize=100))
    cwd = os.getcwd()
    os.chdir(w)
    np.random.seed(seed)
    rt = pd.read_csv("runtime.txt", header=None, index_col=0, sep="\t").drop("output_format").astype(float)
    d = initialize_data(rt)
    enz.append(np.asarray(d["EnzGenes"]).sum(axis=1).mean())
    upt.append(np.asarray(d["UptakeGenes"]).sum(axis=1).mean())
    os.chdir(cwd)
    shutil.rmtree(tmp, ignore_errors=True)
print(f"{os.path.basename(src)}: enzyme genes/taxon mean {np.mean(enz):.3f} (per-seed {np.round(enz, 2).tolist()}); "
      f"transporter genes/taxon mean {np.mean(upt):.3f} (per-seed {np.round(upt, 2).tolist()})")
