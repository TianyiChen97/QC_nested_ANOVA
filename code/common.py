"""Shared config: time points -> DIV folder per run. Outputs go to  nested anova/<timepoint>/ ."""
import os, sys
CODE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(CODE)                       # .../nested anova
BASE = os.path.dirname(ROOT)                       # .../DA data
DIVS = {'early': {'22': 4, '23': 4, '24': 4, '42': 3, '43': 3, '45': 3, 'AR': 4},
        'mid':   {'22': 11, '23': 10, '24': 10, '42': 10, '43': 10, '45': 10, 'AR': 9},
        'late':  {'22': 16, '23': 17, '24': 17, '42': 18, '43': 18, '45': 17, 'AR': 17}}
def setup():
    tp = sys.argv[1]; assert tp in DIVS, f'time point must be one of {list(DIVS)}'
    out = os.path.join(ROOT, tp); os.makedirs(out, exist_ok=True); return tp, out
