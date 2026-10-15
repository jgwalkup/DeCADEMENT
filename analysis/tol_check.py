"""Check the realised drought-tolerance values the model builds for a given input folder (small grid)."""
import os, shutil, sys, tempfile

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
sys.path.insert(0, HERE)
from initialization import initialize_data  # noqa: E402
from sweep_minimal import set_runtime  # noqa: E402

src = os.path.abspath(sys.argv[1])
tmp = tempfile.mkdtemp()
w = os.path.join(tmp, "in")
shutil.copytree(src, w)
set_runtime(os.path.join(w, "runtime.txt"), dict(x=10, y=10, gridsize=100))
os.chdir(w)
np.random.seed(0)
rt = pd.read_csv("runtime.txt", header=None, index_col=0, sep="\t").drop("output_format").astype(float)
d = initialize_data(rt)
tol = np.asarray(d["TaxDroughtTol"], dtype=float)
print(f"{os.path.basename(src)}: n_taxa={int(rt.loc['n_taxa', 1])} unique tolerance values={np.unique(tol.round(6))[:10]} (min {tol.min():.3f}, max {tol.max():.3f}, mean {tol.mean():.3f})")
og = d["OsmoProdConsti"] if "OsmoProdConsti" in d else None
if og is not None:
    s = (np.asarray(d["OsmoProdConsti"]).sum(axis=1) + np.asarray(d["OsmoProdInduci"]).sum(axis=1))
    print("  total osmolyte allocation per taxon (consti+induci): unique", np.unique(s.round(6))[:8])
shutil.rmtree(tmp, ignore_errors=True)
