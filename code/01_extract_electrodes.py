"""Per-electrode mean ISI (<200 ms) + position for every vehicle-control (0uM) well at one time point. Usage: 01_extract_electrodes.py early|mid|late
Active electrode: rate > 0.1 Hz and P90 |neg peak| > 20 uV (same rule as activityscan_metrics.py). 0uM_2 (second plate, Run AR DIV 9) is not used."""
import sys, glob, collections
import h5py, numpy as np, pandas as pd
from common import *
sys.path.insert(0, BASE)
from activityscan_metrics import FR_MIN, AMP_MIN, ISI_MAX
tp, OUT = setup()

def well_electrodes(path, well):
    f = h5py.File(path); R = collections.defaultdict(lambda: dict(T=0., n=0, neg=[], isi=[], xy=None))
    for g in sorted(f['data_store']):
        s = f['data_store'][g]
        if s['well_id'][0] != well: continue
        lsb = s['settings/lsb'][0] * 1e6
        m = s['settings/mapping'][:]; c2m = {c: (e, x, y) for c, e, x, y in m}
        T = (s['stop_time'][0] - s['start_time'][0]) / 1000
        sp = s['spikes'][:]; got = set()
        for ch in np.unique(sp['channel']):
            if ch not in c2m: continue
            x = np.sort(sp[sp['channel'] == ch], order='frameno'); got.add(ch)
            e = c2m[ch]; r = R[e[0]]; r['xy'] = e[1:]
            r['T'] += T; r['n'] += len(x)
            r['neg'].append(-x['amplitude'][x['amplitude'] < 0].astype(float) * lsb)
            r['isi'].append(np.diff(x['frameno']) / 10.0)
        for ch, e in c2m.items():
            if ch not in got:
                r = R[e[0]]; r['xy'] = e[1:]; r['T'] += T
    return R

rows = []
for run, div in DIVS[tp].items():
    path = glob.glob(f"{BASE}/Run_{run}/DIV_{div}/0uM/ActivityScan/*/data.raw.h5")[0]
    print(tp, run, div, path, flush=True)
    for w in range(6):
        R = well_electrodes(path, w); n_rec = len(R)
        for e, r in R.items():
            neg = np.concatenate(r['neg']) if r['neg'] else np.array([]); isi = np.concatenate(r['isi']) if r['isi'] else np.array([])
            fr = r['n'] / r['T']; amp = np.percentile(neg, 90) if len(neg) else 0.
            if fr > FR_MIN and amp > AMP_MIN:
                isi = isi[isi < ISI_MAX]
                if len(isi) == 0: continue
                rows.append(dict(run=run, div=div, well=w + 1, electrode=e, x=r['xy'][0], y=r['xy'][1],
                                 n_spikes=r['n'], rate_hz=fr, n_isi=len(isi), mean_isi_ms=isi.mean(), n_recorded=n_rec))
df = pd.DataFrame(rows); df.to_csv(f'{OUT}/electrodes.csv', index=False)
print(df.groupby(['run', 'well']).size().unstack())
