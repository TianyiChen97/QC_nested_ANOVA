"""All inter-spike intervals (ms, no truncation) of every ACTIVE vehicle-control electrode with >= N_EXTRACT ISIs. Usage: 01_extract_isis.py early|mid|late
Active: rate > 0.1 Hz and P90 |neg peak| > 20 uV (as activityscan_metrics.py). ISIs are computed within each recording config (never across configs), then pooled per electrode."""
import glob, collections, pickle
import h5py, numpy as np, pandas as pd
from common import *
from activityscan_metrics import FR_MIN, AMP_MIN
tp, OUT = setup(); N_EXTRACT = 20

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
            e = c2m[ch]; r = R[e[0]]; r['xy'] = e[1:]; r['T'] += T; r['n'] += len(x)
            r['neg'].append(-x['amplitude'][x['amplitude'] < 0].astype(float) * lsb); r['isi'].append(np.diff(x['frameno']) / 10.0)
        for ch, e in c2m.items():
            if ch not in got: r = R[e[0]]; r['xy'] = e[1:]; r['T'] += T
    return R

meta, isis = [], {}
for run, div in DIVS[tp].items():
    path = glob.glob(f"{BASE}/Run_{run}/DIV_{div}/0uM/ActivityScan/*/data.raw.h5")[0]; print(tp, run, div, flush=True)
    for w in range(6):
        for e, r in well_electrodes(path, w).items():
            neg = np.concatenate(r['neg']) if r['neg'] else np.array([]); isi = np.concatenate(r['isi']) if r['isi'] else np.array([])
            fr = r['n'] / r['T']; amp = np.percentile(neg, 90) if len(neg) else 0.
            if fr > FR_MIN and amp > AMP_MIN and len(isi) >= N_EXTRACT:
                key = f'{run}_w{w+1}_e{e}'; isis[key] = isi.astype(np.float32)
                meta.append(dict(key=key, run=run, well=w + 1, electrode=e, px=r['xy'][0], py=r['xy'][1], n_isi=len(isi), rate_hz=fr))
pd.DataFrame(meta).to_csv(f'{OUT}/electrodes_isi_meta.csv', index=False); pickle.dump(isis, open(f'{OUT}/isis.pkl', 'wb'))
print(len(meta), 'electrodes with >=', N_EXTRACT, 'ISIs')
