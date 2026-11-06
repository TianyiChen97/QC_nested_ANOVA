"""Nested Frechet (2-Wasserstein) ANOVA on per-electrode quantile functions Q_ijk in L2[0,1] (discretised on U).  Usage: 03_frechet_nested_anova.py early|mid|late [N_MIN=30] [gmm|emp]
Barycenter of a set of distributions on R = pointwise mean of their quantile functions. With ||f||^2 = mean_k f(u_k)^2 the decomposition
  SS_total = SS_run + SS_well(run) + SS_e   holds exactly (Hilbert space), same Searle coefficients as the scalar case.
No parametric F reference exists for functional data: tests are permutation tests (run: permute run labels of whole wells; well: permute electrodes among wells within run)."""
import os, sys, pickle
import numpy as np, pandas as pd
from scipy import sparse
from statsmodels.stats.multitest import multipletests
from common import *
tp, OUT0 = setup()
N_MIN = int(sys.argv[2]) if len(sys.argv) > 2 else 30
MODE = sys.argv[3] if len(sys.argv) > 3 else 'gmm'
IN = f'{OUT0}/nmin{N_MIN}'
if os.path.exists(f'{IN}/NOT_ENOUGH_DATA.txt'): print('not enough data'); sys.exit(0)
OUT = IN if MODE == 'gmm' else f'{IN}/sensitivity_empirical_quantiles'; os.makedirs(OUT, exist_ok=True)
NPERM = 5000; rng = np.random.default_rng(1)
th = pd.read_csv(f'{IN}/electrode_fits.csv'); Q = np.load(f'{IN}/Q_gmm.npy' if MODE == 'gmm' else f'{IN}/Q_emp.npy')
th['run'] = th.run.astype(str); th['wid'] = th.run + '_w' + th.well.astype(str)
runs = sorted(th.run.unique(), key=lambda r: (not r.isdigit(), int(r) if r.isdigit() else 0)); wells = sorted(th.wid.unique(), key=lambda w: (runs.index(w.split('_')[0]), w))
ridx = th.run.map({r: i for i, r in enumerate(runs)}).values; widx = th.wid.map({w: i for i, w in enumerate(wells)}).values
well_run = np.array([runs.index(w.split('_')[0]) for w in wells]); N, Kq = Q.shape; A, B = len(runs), len(wells)
sq = lambda D: (D ** 2).mean(-1)                                         # squared L2 norm on [0,1] (grid mean)
Q2 = (Q ** 2).sum(0)

def ss_profiles(widx_, wrun_):
    """returns well means M (B,K), n_ij, run means Mi, n_i, grand G and SS profiles over u"""
    P = sparse.csr_matrix((np.ones(N), (widx_, np.arange(N))), shape=(B, N)); nij = np.asarray(P.sum(1)).ravel(); M = (P @ Q) / nij[:, None]
    ni = np.bincount(wrun_, weights=nij, minlength=A); Mi = np.vstack([(nij[wrun_ == i, None] * M[wrun_ == i]).sum(0) / ni[i] for i in range(A)]); G = Q.mean(0)
    SSA = (ni[:, None] * (Mi - G) ** 2).sum(0); SSB = (nij[:, None] * (M - Mi[wrun_]) ** 2).sum(0); SSE = Q2 - (nij[:, None] * M ** 2).sum(0)
    return dict(M=M, nij=nij, Mi=Mi, ni=ni, G=G, SSA=SSA, SSB=SSB, SSE=SSE)
o = ss_profiles(widx, well_run); nij, ni = o['nij'], o['ni']
dfA, dfB, dfE = A - 1, B - A, N - B
MS = lambda ss, df: ss / df
SSA, SSB, SSE = o['SSA'], o['SSB'], o['SSE']
s1 = (nij ** 2 / ni[well_run]).sum(); k1 = (N - s1) / dfB; k2 = (s1 - (nij ** 2).sum() / N) / dfA; k3 = (N - (ni ** 2).sum() / N) / dfA
def components(SSA, SSB, SSE):
    MSA, MSB, MSE = SSA / dfA, SSB / dfB, SSE / dfE; se = MSE; sb = (MSB - MSE) / k1; sa = (MSA - se - k2 * sb) / k3
    return np.maximum(sa, 0), np.maximum(sb, 0), se
ca, cb, ce = components(SSA, SSB, SSE)                                  # profiles over u
comp = dict(run=ca.mean(), well=cb.mean(), electrode=ce.mean()); tot = sum(comp.values())
FA = (SSA.mean() / dfA) / (SSB.mean() / dfB); FB = (SSB.mean() / dfB) / (SSE.mean() / dfE)

# ---- permutation tests ----
FAp, FBp = [], []; per_stat = np.zeros((NPERM, A)); runs_el = [np.where(ridx == i)[0] for i in range(A)]
SSBi = lambda oo: np.array([(oo['nij'][well_run == i, None] * (oo['M'][well_run == i] - oo['Mi'][i]) ** 2).sum() / Kq for i in range(A)])
SSB_i_obs = SSBi(o); SSE_i_obs = np.array([sq(Q[runs_el[i]] - o['M'][widx[runs_el[i]]]).sum() for i in range(A)])
for b in range(NPERM):
    p = ss_profiles(widx, rng.permutation(well_run)); FAp.append((p['SSA'].mean() / dfA) / (p['SSB'].mean() / dfB))          # run labels of whole wells
    w2 = widx.copy()
    for i in range(A): w2[runs_el[i]] = rng.permutation(widx[runs_el[i]])                                                   # electrodes among wells within run
    p = ss_profiles(w2, well_run); FBp.append((p['SSB'].mean() / dfB) / (p['SSE'].mean() / dfE)); per_stat[b] = SSBi(p)
FAp, FBp = np.array(FAp), np.array(FBp)
pA = (1 + (FAp >= FA).sum()) / (1 + NPERM); pB = (1 + (FBp >= FB).sum()) / (1 + NPERM)
per = []
for i, r in enumerate(runs):
    bi = (well_run == i).sum()
    if bi < 2: per.append(dict(run=r, wells=bi, SS_well=np.nan, eta2=np.nan, p=np.nan)); continue
    eta = SSB_i_obs[i] / (SSB_i_obs[i] + SSE_i_obs[i]); p = (1 + (per_stat[:, i] >= SSB_i_obs[i]).sum()) / (1 + NPERM)
    per.append(dict(run=r, wells=bi, SS_well=SSB_i_obs[i], eta2=eta, p=p))
per = pd.DataFrame(per); ok = per.p.notna(); per.loc[ok, 'p_holm'] = multipletests(per.p[ok], method='holm')[1]

tab = pd.DataFrame([dict(source='Run', SS=SSA.mean(), df=dfA, MS=SSA.mean() / dfA, F_stat=FA, perm_p=pA),
                    dict(source='Well(Run)', SS=SSB.mean(), df=dfB, MS=SSB.mean() / dfB, F_stat=FB, perm_p=pB),
                    dict(source='Electrode (resid)', SS=SSE.mean(), df=dfE, MS=SSE.mean() / dfE, F_stat=np.nan, perm_p=np.nan)])
# descriptive: all eligible electrodes, empirical quantiles, no thinning (variance components only)
Qa = np.load(f'{IN}/Q_emp_all.npy'); ea = pd.read_csv(f'{IN}/eligible_all_meta.csv'); ea['run'] = ea.run.astype(str); ea['wid'] = ea.run + '_w' + ea.well.astype(str); ea = ea[ea.wid.isin(wells)]
Qa = Qa[ea.index.values] if len(Qa) == len(pd.read_csv(f'{IN}/eligible_all_meta.csv')) else Qa
def comps_generic(Qm, run_i, well_i):
    n_ = len(Qm); Bm = well_i.max() + 1; Ar = run_i.max() + 1
    P = sparse.csr_matrix((np.ones(n_), (well_i, np.arange(n_))), shape=(Bm, n_)); nj = np.asarray(P.sum(1)).ravel(); Mm = (P @ Qm) / nj[:, None]
    wr = np.zeros(Bm, int); wr[well_i] = run_i; nI = np.bincount(wr, weights=nj, minlength=Ar); MI = np.vstack([(nj[wr == i, None] * Mm[wr == i]).sum(0) / nI[i] for i in range(Ar)]); Gm = Qm.mean(0)
    sA = (nI[:, None] * (MI - Gm) ** 2).sum(0).mean(); sB = (nj[:, None] * (Mm - MI[wr]) ** 2).sum(0).mean(); sE = ((Qm ** 2).sum(0) - (nj[:, None] * Mm ** 2).sum(0)).mean()
    d1, d2, d3 = Ar - 1, Bm - Ar, n_ - Bm; ss1 = (nj ** 2 / nI[wr]).sum(); K1 = (n_ - ss1) / d2; K2 = (ss1 - (nj ** 2).sum() / n_) / d1; K3 = (n_ - (nI ** 2).sum() / n_) / d1
    se = sE / d3; sb = (sB / d2 - se) / K1; sa = (sA / d1 - se - K2 * sb) / K3; return dict(run=max(sa, 0), well=max(sb, 0), electrode=se)
ea_run = ea.run.map({r: i for i, r in enumerate(runs)}).values; ea_w = ea.wid.map({w: i for i, w in enumerate(wells)}).values
desc = comps_generic(Qa, ea_run, ea_w)

L = []; add = L.append; sz = pd.Series(widx).value_counts()
add(f'# Nested Frechet (W2) ANOVA of per-electrode ISI distributions — {tp} DIV {sorted(set(DIVS[tp].values()))}, vehicle control\n')
add(f"Per-electrode distribution: " + ('2-component GMM on log10(ISI ms), converted to quantile function' if MODE == 'gmm' else 'EMPIRICAL quantile function of log10(ISI ms) (sensitivity: no mixture model)') +
    f"; electrodes need >= {N_MIN} ISIs; thinned to >= {D_BY_TP[tp]:.0f} um; wells kept with >= 10 thinned electrodes: {B} wells in {A} runs, N = {N} electrodes (m per well {sz.min()}-{sz.max()}). Quantile grid K = {Kq}.\n")
add('## Frechet ANOVA table (SS in (log10 ISI)^2 units, L2 on quantile functions)\n'); add(tab.round(5).to_markdown(index=False))
add(f'\nF_stat is the usual ratio MS_run/MS_well(run) and MS_well(run)/MS_e; its null distribution is NOT F, p-values are from {NPERM} permutations.\n')
add('## Per-run well test (permutation; eta2 = SS_well / (SS_well + SS_e) within run, Holm across runs)\n'); add(per.round(5).to_markdown(index=False))
add('\n## Variance components (Frechet variance, (log10 ISI)^2) and share of total\n')
vc = pd.DataFrame([dict(method='thinned (tested)', **{f'var_{k}': v for k, v in comp.items()}, **{f'prop_{k}': v / tot for k, v in comp.items()}, N=N),
                   dict(method='all eligible electrodes, empirical Q (descriptive)', **{f'var_{k}': v for k, v in desc.items()}, **{f'prop_{k}': v / sum(desc.values()) for k, v in desc.items()}, N=len(Qa))])
add(vc.round(5).to_markdown(index=False)); vc.to_csv(f'{OUT}/variance_components.csv', index=False); tab.to_csv(f'{OUT}/frechet_anova_table.csv', index=False); per.to_csv(f'{OUT}/per_run_well_tests.csv', index=False)
pd.DataFrame(dict(u=U, var_run=ca, var_well=cb, var_electrode=ce)).to_csv(f'{OUT}/variance_profile_by_quantile_level.csv', index=False)
open(f'{OUT}/results.md', 'w').write('\n'.join(L)); print('\n'.join(L))
pickle.dump(dict(th=th, Q=Q, runs=runs, wells=wells, o=o, tab=tab, per=per, comp=comp, ca=ca, cb=cb, ce=ce, pA=pA, pB=pB, FA=FA, FB=FB, vc=vc, MODE=MODE, N_MIN=N_MIN, k=(k1, k2, k3), tp=tp, widx=widx, well_run=well_run), open(f'{OUT}/_res.pkl', 'wb'))
