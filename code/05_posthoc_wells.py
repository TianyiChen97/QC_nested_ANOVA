"""Post-hoc: which runs have wells whose means differ within run, and which well pairs. y = ln(mean ISI), thinned data. Usage: 05_posthoc_wells.py early|mid|late
(1) per-run omnibus F with pooled MSE_e, Holm across runs. (2) Tukey-Kramer pairwise wells within run (pooled MSE_e). Sensitivity: run-specific MSE one-way ANOVA."""
import pickle, itertools
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from common import *
tp, OUT = setup()
R = pickle.load(open(f'{OUT}/_res.pkl', 'rb')); th, res = R['th'].copy(), R['res']
th['y'] = np.log(th.mean_isi_ms); MSE, dfE = res['MSE'], res['dfE']
runs = sorted(th.run.unique(), key=lambda r: (not r.isdigit(), int(r) if r.isdigit() else 0))
om, pw = [], []
for r in runs:
    g = th[th.run == r]; s = g.groupby('well').y.agg(['mean', 'size']); k = len(s)
    if k < 2: om.append(dict(run=r, wells=k, F=np.nan, p=np.nan)); continue
    ss = float((s['size'] * (s['mean'] - g.y.mean()) ** 2).sum()); F = ss / (k - 1) / MSE
    om.append(dict(run=r, wells=k, F=F, p=stats.f.sf(F, k - 1, dfE)))
    for a, b in itertools.combinations(s.index, 2):
        d = s.loc[a, 'mean'] - s.loc[b, 'mean']; se = np.sqrt(MSE / 2 * (1 / s.loc[a, 'size'] + 1 / s.loc[b, 'size']))
        q = abs(d) / se; p = stats.studentized_range.sf(q, k, dfE); h = stats.studentized_range.ppf(.95, k, dfE) * se
        pw.append(dict(run=r, well_a=a, well_b=b, n_a=int(s.loc[a, 'size']), n_b=int(s.loc[b, 'size']), diff_log=d, ci_lo=d - h, ci_hi=d + h, ratio=np.exp(d), p_tukey=p))
om = pd.DataFrame(om); ok = om.p.notna(); om.loc[ok, 'p_holm'] = multipletests(om.p[ok], method='holm')[1]; om['sig_holm'] = om.p_holm < .05
pw = pd.DataFrame(pw); pw['sig'] = pw.p_tukey < .05
sens = []
for r in runs:
    g = th[th.run == r]
    if g.well.nunique() < 2: continue
    f = stats.f_oneway(*[x.y.values for _, x in g.groupby('well') if len(x) > 1]); sens.append(dict(run=r, F_own_MSE=f.statistic, p_own_MSE=f.pvalue))
om = om.merge(pd.DataFrame(sens), on='run', how='left')
om.to_csv(f'{OUT}/posthoc_run_omnibus.csv', index=False); pw.to_csv(f'{OUT}/posthoc_well_pairs.csv', index=False)
pd.set_option('display.width', 200); print(om.round(4).to_string(index=False)); print(pw[pw.sig].round(4).to_string(index=False))
fig, axs = plt.subplots(1, 7, figsize=(18, 5.2), sharey=True); tab = plt.get_cmap('tab10')
for i, (ax, r) in enumerate(zip(axs, runs)):
    g = th[th.run == r]; s = g.groupby('well').y.agg(['mean', 'size']); ci = stats.t.ppf(.975, dfE) * np.sqrt(MSE / s['size'])
    ax.errorbar(s.index, s['mean'], ci, fmt='o', color=tab(i), capsize=3); ax.axhline(g.y.mean(), color='grey', ls=':', lw=.8)
    ax.set_xticks(range(1, 7)); ax.set_xlim(.4, 6.6); ax.set_xlabel('well')
    o = om[om.run == r].iloc[0]; sg = pw[(pw.run == r) & pw.sig].sort_values('p_tukey')
    t = f"{int(o.wells)} wells\n" + ('single well' if np.isnan(o.F) else f"F={o.F:.2f}, Holm p={o.p_holm:.3g}" + ('  *' if o.sig_holm else ''))
    t += ''.join(f"\nw{a.well_a}≠w{a.well_b} p={a.p_tukey:.3g}" for a in sg.head(7).itertuples()) + (f"\n+{len(sg)-7} more" if len(sg) > 7 else '')
    ax.set_title(f'run {r}', color=tab(i), fontweight='bold'); ax.text(.5, -.3, t, transform=ax.transAxes, ha='center', va='top', fontsize=7.5)
axs[0].set_ylabel('well mean  ln(mean ISI)  ±95% CI'); fig.suptitle(f'{tp}: well means within each run (Tukey–Kramer pairs, pooled MS_e)', y=.99); fig.subplots_adjust(bottom=.45)
fig.savefig(f'{OUT}/posthoc_wells_within_run.png', dpi=160, bbox_inches='tight', facecolor='white')
