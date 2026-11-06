"""Nested random-effects ANOVA run > well > electrode on ln(mean ISI<200 ms), vehicle control, one time point.
Usage: 03_nested_anova.py early|mid|late [D_MIN_um=100] [MIN_M=10] [log|raw]
Spatial thinning at D_MIN um (electrode position only, response never used). Unbalanced -> Searle coefficients + REML cross-check
+ permutation tests (parametric F-test and permutation test are both reported)."""
import os, sys, pickle
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from common import *
tp, OUT = setup()
D_MIN = float(sys.argv[2]) if len(sys.argv) > 2 else 100.0
MIN_M = int(sys.argv[3]) if len(sys.argv) > 3 else 10
TRANSFORM = sys.argv[4] if len(sys.argv) > 4 else 'log'
if TRANSFORM == 'raw': OUT = f'{OUT}/sensitivity_raw'; os.makedirs(OUT, exist_ok=True)
NPERM, SEED = 5000, 0
rng = np.random.default_rng(SEED)
IN = f'{ROOT}/{tp}'

df = pd.read_csv(f'{IN}/electrodes.csv').rename(columns={'x': 'px', 'y': 'py'})   # px,py = position (um); y = response
df['y'] = np.log(df.mean_isi_ms) if TRANSFORM == 'log' else df.mean_isi_ms
df['run'] = df.run.astype(str); df['wid'] = df.run + '_w' + df.well.astype(str)

def thin(g):
    """greedy: random visiting order, keep an electrode if no kept one is within D_MIN (positions only)"""
    xy = g[['px', 'py']].values; kept = []
    for i in rng.permutation(len(g)):
        if all(np.hypot(*(xy[i] - xy[k])) >= D_MIN for k in kept): kept.append(i)
    return g.iloc[sorted(kept)]
thinned = pd.concat([thin(g) for _, g in df.groupby('wid')])
cnt = pd.DataFrame({'active_all': df.groupby('wid').size(), 'm_thinned': thinned.groupby('wid').size()}).fillna(0).astype(int)
keep = cnt.index[cnt.m_thinned >= MIN_M]
print(cnt.sort_values('m_thinned').T.to_string()); print('wells kept:', len(keep), 'of', len(cnt))
th = thinned[thinned.wid.isin(keep)].copy()
th.to_csv(f'{OUT}/electrodes_thinned.csv', index=False); cnt.to_csv(f'{OUT}/well_counts.csv')

def anova(d, label):
    d = d.copy(); N = len(d)
    nij = d.groupby(['run', 'wid']).size(); ni = d.groupby('run').size()
    a, B = d.run.nunique(), d.wid.nunique()
    gm = d.y.mean(); mi = d.groupby('run').y.mean(); mij = d.groupby(['run', 'wid']).y.mean()
    SSA = float((ni * (mi - gm) ** 2).sum())
    SSB = float((nij * (mij - mi.reindex(mij.index.get_level_values(0)).values) ** 2).sum())
    SSE = float(((d.y - d.groupby('wid').y.transform('mean')) ** 2).sum())
    dfA, dfB, dfE = a - 1, B - a, N - B
    MSA, MSB, MSE = SSA / dfA, SSB / dfB, SSE / dfE
    s1 = (nij ** 2 / ni.reindex(nij.index.get_level_values(0)).values).sum()
    k1 = (N - s1) / dfB; k2 = (s1 - (nij ** 2).sum() / N) / dfA; k3 = (N - (ni ** 2).sum() / N) / dfA
    s2e = MSE; s2b = (MSB - MSE) / k1; s2a = (MSA - s2e - k2 * s2b) / k3          # Searle: E[MSA]=s2e+k2*s2b+k3*s2a ; E[MSB]=s2e+k1*s2b
    comp = {'run': max(s2a, 0), 'well': max(s2b, 0), 'electrode': s2e}
    FB = MSB / MSE; pB = stats.f.sf(FB, dfB, dfE); FA = MSA / MSB; pA = stats.f.sf(FA, dfA, dfB)
    per = []
    for r, g in d.groupby('run'):
        nj = g.groupby('wid').size(); mj = g.groupby('wid').y.mean(); b = len(nj)
        if b < 2: per.append(dict(run=r, wells=b, F=np.nan, df1=b - 1, p=np.nan)); continue
        ss = float((nj * (mj - g.y.mean()) ** 2).sum()); F = ss / (b - 1) / MSE
        per.append(dict(run=r, wells=b, F=F, df1=b - 1, p=stats.f.sf(F, b - 1, dfE)))
    tab = pd.DataFrame([dict(source='Run', SS=SSA, df=dfA, MS=MSA, F=FA, p=pA, den='MS_well(run)'),
                        dict(source='Well(Run)', SS=SSB, df=dfB, MS=MSB, F=FB, p=pB, den='MS_electrode'),
                        dict(source='Electrode (resid)', SS=SSE, df=dfE, MS=MSE, F=np.nan, p=np.nan, den='')])
    tot = sum(comp.values())
    return dict(label=label, N=N, a=a, B=B, comp=comp, prop={k: v / tot for k, v in comp.items()}, tab=tab, per=pd.DataFrame(per),
                MSA=MSA, MSB=MSB, MSE=MSE, dfA=dfA, dfB=dfB, dfE=dfE, k=(k1, k2, k3))

def reml(d):
    m = smf.mixedlm('y ~ 1', d, groups='run', re_formula='1', vc_formula={'well': '0 + C(wid)'}).fit(reml=True)
    c = {'run': float(m.cov_re.iloc[0, 0]), 'well': float(m.vcomp[0]), 'electrode': float(m.scale)}; tot = sum(c.values())
    return dict(comp=c, prop={k: v / tot for k, v in c.items()}, mu=float(m.fe_params.iloc[0]))

def perm_tests(d, obs):
    """run test: permute run labels across wells (well sizes kept). well test: permute electrodes among wells within run."""
    wells = d.groupby('wid').agg(run=('run', 'first')).reset_index(); base = d[['run', 'wid', 'y']].copy(); FA, FB = [], []
    for _ in range(NPERM):
        w2 = wells.copy(); w2['run'] = rng.permutation(w2.run.values)
        dd = base.drop(columns='run').merge(w2[['wid', 'run']], on='wid')
        try: FA.append(anova(dd, '')['tab'].F.iloc[0])
        except Exception: pass
        dd = base.copy(); dd['wid'] = dd.groupby('run').wid.transform(lambda s: rng.permutation(s.values))
        try: FB.append(anova(dd, '')['tab'].F.iloc[1])
        except Exception: pass
    FA, FB = np.array(FA), np.array(FB)
    return dict(pA=(1 + (FA >= obs['tab'].F.iloc[0]).sum()) / (1 + len(FA)), pB=(1 + (FB >= obs['tab'].F.iloc[1]).sum()) / (1 + len(FB)))

allk = df[df.wid.isin(keep)]
res_th = anova(th, 'thinned'); res_all = anova(allk, 'all electrodes (same wells)')
rm_th = reml(th); rm_all = reml(allk); pt = perm_tests(th, res_th)

P = []; add = P.append; sz = th.groupby('wid').size()
unit = 'ln(mean ISI [ms])' if TRANSFORM == 'log' else 'mean ISI [ms] (raw)'
add(f'# Nested ANOVA: run > well > electrode, {tp} DIV {sorted(set(DIVS[tp].values()))}, vehicle control, y = {unit}\n')
add(f'Thinning: min pairwise distance d = {D_MIN:.0f} um; wells kept with m >= {MIN_M} thinned electrodes: {len(keep)} of {len(cnt)} wells, '
    f'{th.run.nunique()} runs, N = {len(th)} electrodes. Unbalanced design (m per well {sz.min()}-{sz.max()}).\n')
add('## ANOVA table (thinned); p = parametric F-test\n'); add(res_th['tab'].round(4).to_markdown(index=False))
add(f"\nPermutation p ({NPERM} perms): run = {pt['pA']:.4f}, well(run) = {pt['pB']:.4f}\n")
add('## Per-run well test (thinned; parametric F with pooled MS_electrode)\n'); add(res_th['per'].round(4).to_markdown(index=False))
add('\n## Variance components\n'); rows = []
for nm, r, n in [('ANOVA thinned', res_th, len(th)), ('ANOVA all electrodes', res_all, len(allk)), ('REML thinned', rm_th, len(th)), ('REML all electrodes', rm_all, len(allk))]:
    rows.append(dict(method=nm, **{f'var_{k}': v for k, v in r['comp'].items()}, **{f'prop_{k}': v for k, v in r['prop'].items()}, N=n))
vc = pd.DataFrame(rows); add(vc.round(4).to_markdown(index=False)); vc.to_csv(f'{OUT}/variance_components.csv', index=False)
res_th['tab'].to_csv(f'{OUT}/anova_table_thinned.csv', index=False)
open(f'{OUT}/results.md', 'w').write('\n'.join(P)); print('\n'.join(P))
pickle.dump(dict(th=th, res=res_th, vc=vc, pt=pt, cnt=cnt, keep=list(keep), tp=tp, D_MIN=D_MIN, MIN_M=MIN_M, transform=TRANSFORM), open(f'{OUT}/_res.pkl', 'wb'))
