"""Matched-electrode comparison of the scalar mean-ISI nested ANOVA and the Frechet (quantile-function) nested ANOVA.   Usage: 05_matched_checks.py early|mid|late [N_MIN=30]
Electrode set S = the thinned Frechet electrodes (>= N_MIN ISIs, nmin<N>/electrode_fits.csv) that also have >= 10 ISIs < 200 ms. Same electrodes in all four analyses:
  A  scalar  y = ln(mean ISI | ISI<200 ms)      (the original mean-ISI definition)
  A' scalar  y = ln(mean of ALL ISIs)           (no truncation)
  B  Frechet quantile function of log10 ISI, ISIs < 200 ms only (empirical quantiles)
  B' Frechet quantile function of log10 ISI, all ISIs (empirical quantiles)
All tests are permutation tests (5000); for the scalar rows the parametric F-test p is also given."""
import pickle
import numpy as np, pandas as pd
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests
from common import *
tp, OUT0 = setup(); N_MIN = int(sys.argv[2]) if len(sys.argv) > 2 else 30; IN = f'{OUT0}/nmin{N_MIN}'; OUT = f'{IN}/matched_checks'; os.makedirs(OUT, exist_ok=True)
NPERM = 5000; isis = pickle.load(open(f'{OUT0}/isis.pkl', 'rb'))
th = pd.read_csv(f'{IN}/electrode_fits.csv'); th['run'] = th.run.astype(str); th['wid'] = th.run + '_w' + th.well.astype(str)
tr = {k: isis[k].astype(float) for k in th.key}; n200 = np.array([(tr[k] < 200).sum() for k in th.key]); S = n200 >= 10
th = th[S].reset_index(drop=True); print(f'{tp}: {len(th)} of {S.size} thinned electrodes have >=10 ISIs<200ms; dropped {int((~S).sum())}')
cnt = th.wid.value_counts(); th = th[th.wid.map(cnt) >= 10].reset_index(drop=True)          # re-apply the >=10 electrodes per well rule
runs = sorted(th.run.unique(), key=lambda r: (not r.isdigit(), int(r) if r.isdigit() else 0)); wells = sorted(th.wid.unique(), key=lambda w: (runs.index(w.split('_')[0]), w))
ridx = th.run.map({r: i for i, r in enumerate(runs)}).values; widx = th.wid.map({w: i for i, w in enumerate(wells)}).values
well_run = np.array([runs.index(w.split('_')[0]) for w in wells]); N = len(th); A, B = len(runs), len(wells); print(f'  analysed: N={N}, wells={B}, runs={A}')
lg = lambda x: np.log10(np.clip(x, 0.1, None))
y_trunc = np.log([tr[k][tr[k] < 200].mean() for k in th.key])[:, None]; y_full = np.log([tr[k].mean() for k in th.key])[:, None]
Qt = np.vstack([np.quantile(lg(tr[k][tr[k] < 200]), U) for k in th.key]); Qf = np.vstack([np.quantile(lg(tr[k]), U) for k in th.key])
runs_el = [np.where(ridx == i)[0] for i in range(A)]
def analyze(Q, seed):
    rng = np.random.default_rng(seed); d = Q.shape[1]; Q2 = (Q ** 2).sum(0)
    def prof(w_, wr_):
        P = sparse.csr_matrix((np.ones(N), (w_, np.arange(N))), shape=(B, N)); nij = np.asarray(P.sum(1)).ravel(); M = (P @ Q) / nij[:, None]
        ni = np.bincount(wr_, weights=nij, minlength=A); Mi = np.vstack([(nij[wr_ == i, None] * M[wr_ == i]).sum(0) / ni[i] for i in range(A)]); G = Q.mean(0)
        return dict(M=M, nij=nij, ni=ni, Mi=Mi, SSA=((ni[:, None] * (Mi - G) ** 2).sum(0)).mean(), SSB=((nij[:, None] * (M - Mi[wr_]) ** 2).sum(0)).mean(), SSE=(Q2 - (nij[:, None] * M ** 2).sum(0)).mean())
    o = prof(widx, well_run); nij, ni = o['nij'], o['ni']; dfA, dfB, dfE = A - 1, B - A, N - B
    s1 = (nij ** 2 / ni[well_run]).sum(); k1 = (N - s1) / dfB; k2 = (s1 - (nij ** 2).sum() / N) / dfA; k3 = (N - (ni ** 2).sum() / N) / dfA
    MSA, MSB, MSE = o['SSA'] / dfA, o['SSB'] / dfB, o['SSE'] / dfE; sb = (MSB - MSE) / k1; sa = (MSA - MSE - k2 * sb) / k3
    comp = dict(run=max(sa, 0), well=max(sb, 0), electrode=MSE); tot = sum(comp.values()); FA, FB = MSA / MSB, MSB / MSE
    fa, fb = [], []; per_stat = np.zeros((NPERM, A))
    SSBi = lambda oo: np.array([(oo['nij'][well_run == i, None] * (oo['M'][well_run == i] - oo['Mi'][i]) ** 2).sum() / d for i in range(A)])
    for b in range(NPERM):
        p = prof(widx, rng.permutation(well_run)); fa.append((p['SSA'] / dfA) / (p['SSB'] / dfB))
        w2 = widx.copy()
        for i in range(A): w2[runs_el[i]] = rng.permutation(widx[runs_el[i]])
        p = prof(w2, well_run); fb.append((p['SSB'] / dfB) / (p['SSE'] / dfE)); per_stat[b] = SSBi(p)
    pA = (1 + (np.array(fa) >= FA).sum()) / (1 + NPERM); pB = (1 + (np.array(fb) >= FB).sum()) / (1 + NPERM)
    obs = SSBi(o); pr = np.array([np.nan if (well_run == i).sum() < 2 else (1 + (per_stat[:, i] >= obs[i]).sum()) / (1 + NPERM) for i in range(A)]); ok = ~np.isnan(pr)
    ph = np.full(A, np.nan); ph[ok] = multipletests(pr[ok], method='holm')[1]
    out = dict(pA=pA, pB=pB, pr=comp['run'] / tot, pw=comp['well'] / tot, pe=comp['electrode'] / tot, flagged=','.join(r for r, p in zip(runs, ph) if p < .05))
    if d == 1: out.update(F_pA=stats.f.sf(FA, dfA, dfB), F_pB=stats.f.sf(FB, dfB, dfE))
    return out
rows = []
for name, Q, seed in [("A   scalar ln mean ISI (<200 ms)", y_trunc, 1), ("A'  scalar ln mean ISI (all ISIs)", y_full, 2),
                      ("B   Frechet, ISIs<200 ms only", Qt, 3), ("B'  Frechet, all ISIs (empirical Q)", Qf, 4)]:
    r = analyze(Q, seed); r['analysis'] = name; rows.append(r); print(name, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
df = pd.DataFrame(rows).set_index('analysis'); df.insert(0, 'N', N); df.insert(1, 'wells', B); df.to_csv(f'{OUT}/matched_checks.csv'); open(f'{OUT}/matched_checks.md', 'w').write(df.round(4).to_markdown())
