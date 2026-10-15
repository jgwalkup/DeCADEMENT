"""Community-mean per-cell production totals (sum over genes) that a given input folder produces after initialisation."""
import os, shutil, sys, tempfile

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
sys.path.insert(0, HERE)
from initialization import initialize_data  # noqa: E402
from sweep_minimal import set_runtime  # noqa: E402

src = os.path.abspath(sys.argv[1])
side = int(sys.argv[2]) if len(sys.argv) > 2 else 10
tmp = tempfile.mkdtemp()
w = os.path.join(tmp, "in")
shutil.copytree(src, w)
set_runtime(os.path.join(w, "runtime.txt"), dict(x=side, y=side, gridsize=side * side))
os.chdir(w)
np.random.seed(int(sys.argv[3]) if len(sys.argv) > 3 else 0)
rt = pd.read_csv("runtime.txt", header=None, index_col=0, sep="\t").drop("output_format").astype(float)
d = initialize_data(rt)
print(os.path.basename(src), "keys:", [k for k in d if any(s in k for s in ("Enz", "Uptake", "HSP", "Osmo"))])
def rowsum(name):
    v = d[name]
    return np.asarray(v.sum(axis=1) if hasattr(v, "sum") and getattr(v, "ndim", 1) > 1 else v, dtype=float)
for name in ("EnzProdConstit", "EnzProdInduce", "UptakeGenesCost", "HSPProdConsti", "HSPProdInduci", "OsmoProdConsti", "OsmoProdInduci"):
    if name in d:
        s = rowsum(name)
        print(f"  {name:18s} per-cell total: mean {s.mean():.6f}  sd {s.std():.6f}  min {s.min():.6f}  max {s.max():.6f}")
shutil.rmtree(tmp, ignore_errors=True)
