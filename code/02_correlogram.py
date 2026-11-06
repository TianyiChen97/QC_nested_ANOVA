"""Residual spatial correlogram of ln(mean ISI) (residual = value - well mean), pooled over wells. Usage: 02_correlogram.py early|mid|late"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from common import *
tp, OUT = setup()
df = pd.read_csv(f'{OUT}/electrodes.csv'); df['ly'] = np.log(df.mean_isi_ms)
df['res'] = df.ly - df.groupby(['run', 'well']).ly.transform('mean')
bins = [0, 20, 40, 60, 100, 150, 250, 400, 700, 1000, 1500]
num = np.zeros(len(bins) - 1); cnt = np.zeros(len(bins) - 1)
for _, g in df.groupby(['run', 'well']):
    if len(g) < 5: continue
    e = g.res.values; xy = g[['x', 'y']].values; i, j = np.triu_indices(len(g), 1)
    d = np.hypot(*(xy[i] - xy[j]).T); p = e[i] * e[j]
    for b in range(len(bins) - 1):
        m = (d >= bins[b]) & (d < bins[b + 1]); num[b] += p[m].sum(); cnt[b] += m.sum()
v = (df.res ** 2).mean(); rho = num / np.maximum(cnt, 1) / v
out = pd.DataFrame(dict(d_lo=bins[:-1], d_hi=bins[1:], pairs=cnt.astype(int), rho=rho.round(4))); out.to_csv(f'{OUT}/correlogram.csv', index=False); print(out.to_string(index=False))
fig, ax = plt.subplots(figsize=(5, 3.2)); ax.plot(np.add(bins[:-1], bins[1:]) / 2, rho, 'o-'); ax.axhline(0, color='grey', lw=.7)
ax.axvline(100, color='r', ls=':', lw=.8); ax.set_xscale('log'); ax.set_xlabel('pair distance (µm)'); ax.set_ylabel('residual correlation'); ax.set_title(f'{tp}: residual correlogram (red = d used)')
fig.savefig(f'{OUT}/correlogram.png', dpi=150, bbox_inches='tight')
