#!/usr/bin/env python3
"""Run perturbation_parallel_template.ipynb non-interactively: executes its code cells top to bottom in one namespace, from this folder.
(nbconvert / papermill are not installed in the scenicplus environment; the notebook has no magics, so plain exec is equivalent.)

    PERTURB_INPUT_DIR=<results dir> PERTURB_CONFIGS=config3 PERTURB_H5MU=scplusmdata_perturb.h5mu PERTURB_N_JOBS=40 \
        python run_notebook.py [notebook.ipynb]
"""
import json
import os
import sys
import time

nb_path = sys.argv[1] if len(sys.argv) > 1 else "perturbation_parallel_template.ipynb"
os.chdir(os.path.dirname(os.path.abspath(__file__)))
cells = [c for c in json.load(open(nb_path))["cells"] if c["cell_type"] == "code"]
ns = {"__name__": "__main__"}
for i, c in enumerate(cells):
    src = "".join(c["source"])
    assert not any(l.lstrip().startswith(("%", "!")) for l in src.splitlines()), f"cell {i} has a magic/shell line"
    t0 = time.time()
    exec(compile(src, f"{nb_path}[code cell {i}]", "exec"), ns)
    print(f"--- code cell {i} done ({time.time() - t0:.0f}s)", flush=True)
