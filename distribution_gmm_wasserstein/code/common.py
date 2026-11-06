"""Shared config. Reuses time-point -> DIV map from ../../code/common.py. Outputs go to  distribution_gmm_wasserstein/<timepoint>/ ."""
import os, sys
import numpy as np
CODE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(CODE)                      # .../distribution_gmm_wasserstein
NA = os.path.dirname(ROOT)                        # .../nested anova
BASE = os.path.dirname(NA)                        # .../DA data
sys.path.insert(0, BASE)                          # for activityscan_metrics
import importlib.util
_s = importlib.util.spec_from_file_location('na_common', os.path.join(NA, 'code', 'common.py')); _m = importlib.util.module_from_spec(_s); _s.loader.exec_module(_m)
DIVS = _m.DIVS                                    # time point -> DIV folder per run (nested anova/code/common.py)
D_BY_TP = {'early': 100.0, 'mid': 150.0, 'late': 150.0}   # um, from the residual correlograms in nested anova/<tp>/correlogram.csv
K = 200                                           # quantile grid size: u_k = (k-0.5)/K
U = (np.arange(1, K + 1) - 0.5) / K
def setup():
    tp = sys.argv[1]; assert tp in DIVS, f'time point must be one of {list(DIVS)}'
    out = os.path.join(ROOT, tp); os.makedirs(out, exist_ok=True); return tp, out
