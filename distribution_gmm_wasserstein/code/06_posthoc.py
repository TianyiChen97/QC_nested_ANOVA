"""Post hoc for the nested Frechet ANOVA: which wells / runs differ.   Usage: 06_posthoc.py early|mid|late [N_MIN=30] [gmm|emp]
WELLS within a run (family = that run; null: wells of the run are exchangeable, so electrodes are permuted among the run's wells; same scan protocol inside a run):
   pairwise  T_ab = n_a n_b/(n_a+n_b) * ||Qbar_a - Qbar_b||^2 / MS_e        (n-weighted W2^2 between well barycenters, studentised like q^2)
   one-vs-rest T_a = n_a (n_r-n_a)/n_r * ||Qbar_a - Qbar_(run minus a)||^2 / MS_e
   adjusted p = max-T permutation (Westfall-Young) over the family, 5000 perms.
RUNS (family = all runs; units = wells; null: wells exchangeable across runs, so run labels of whole wells are permuted):
   pairwise  T_ij = b_i b_j/(b_i+b_j) * ||m_i - m_j||^2 / MS_well(run),  m_i = UNweighted mean of the run's well barycenters, b_i = number of wells
   one-vs-rest T_i = b_i (B-b_i)/B * ||m_i - m_(others)||^2 / MS_well(run).   Runs with 1 well have no within-run replication (their T is driven by a single well)."""
import pickle, itertools
import numpy as np, pandas as pd
from scipy import sparse
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from common import *
tp, OUT0 = setup(); N_MIN = int(sys.argv[2]) if len(sys.argv) > 2 else 30; MODE = sys.argv[3] if len(sys.argv) > 3 else 'gmm'
OUT = f'{OUT0}/nmin{N_MIN}' + ('' if MODE == 'gmm' else '/sensitivity_empirical_quantiles'); PH = f'{OUT}/posthoc'; os.makedirs(PH, exist_ok=True)
R = pickle.load(open(f'{OUT}/_res.pkl', 'rb')); Q, runs, wells, o, widx, well_run, tab = R['Q'], R['runs'], R['wells'], R['o'], R['widx'], R['well_run'], R['tab'].set_index('source')
A, B, (N, K) = len(runs), len(wells), Q.shape; MSE = tab.loc['Electrode (resid)', 'MS']; MSB = tab.loc['Well(Run)', 'MS']; NP = 5000; rng = np.random.default_rng(7)
sq = lambda D: (D ** 2).mean(-1); nij = o['nij']; M = o['M']; wl = [w.split('_w')[1] for w in wells]
runs_el = [np.where(well_run[widx] == i)[0] for i in range(A)]; run_wells = [np.where(well_run == i)[0] for i in range(A)]
nrun = np.array([nij[w].sum() for w in run_wells])
# ---------- wells within run ----------
pa, pb, prun = [], [], []
for i in range(A):
    for a, b in itertools.combinations(run_wells[i], 2): pa.append(a); pb.append(b); prun.append(i)
pa, pb, prun = map(np.array, (pa, pb, prun)); cpair = nij[pa] * nij[pb] / (nij[pa] + nij[pb])
rest_ok = np.array([len(run_wells[well_run[j]]) >= 2 for j in range(B)]); wrun = well_run; crest = nij * (nrun[wrun] - nij) / nrun[wrun]
def well_stats(Mw):
    Tp = cpair * sq(Mw[pa] - Mw[pb]) / MSE
    S = np.vstack([(nij[run_wells[i], None] * Mw[run_wells[i]]).sum(0) for i in range(A)]); den = np.where(rest_ok, nrun[wrun] - nij, 1)
    rest = (S[wrun] - nij[:, None] * Mw) / den[:, None]; Tr = np.where(rest_ok, crest * sq(Mw - rest) / MSE, 0.0); return Tp, Tr
Tp_obs, Tr_obs = well_stats(M); mp = np.zeros((NP, A)); mr = np.zeros((NP, A))
for b in range(NP):
    w2 = widx.copy()
    for i in range(A): w2[runs_el[i]] = rng.permutation(widx[runs_el[i]])
    P = sparse.csr_matrix((np.ones(N), (w2, np.arange(N))), shape=(B, N)); Mw = (P @ Q) / nij[:, None]; Tp, Tr = well_stats(Mw)
    for i in range(A):
        sel = prun == i
        if sel.any(): mp[b, i] = Tp[sel].max()
        if len(run_wells[i]) >= 2: mr[b, i] = Tr[run_wells[i]].max()
padj_pair = np.array([(1 + (mp[:, prun[k]] >= Tp_obs[k]).sum()) / (1 + NP) for k in range(len(pa))])
padj_rest = np.array([np.nan if not rest_ok[j] else (1 + (mr[:, wrun[j]] >= Tr_obs[j]).sum()) / (1 + NP) for j in range(B)])
rest_means = None
wp = pd.DataFrame(dict(run=[runs[i] for i in prun], well_a=[wl[a] for a in pa], well_b=[wl[b] for b in pb], n_a=nij[pa].astype(int), n_b=nij[pb].astype(int),
                       W2_distance=np.sqrt(sq(M[pa] - M[pb])), T=Tp_obs, p_adj=padj_pair)); wp['sig'] = wp.p_adj < .05
Srun = np.vstack([(nij[run_wells[i], None] * M[run_wells[i]]).sum(0) for i in range(A)]); restM = (Srun[wrun] - nij[:, None] * M) / np.where(rest_ok, nrun[wrun] - nij, 1)[:, None]
wr = pd.DataFrame(dict(run=[runs[i] for i in wrun], well=wl, n=nij.astype(int), W2_to_rest_of_run=np.where(rest_ok, np.sqrt(sq(M - restM)), np.nan), T=np.where(rest_ok, Tr_obs, np.nan), p_adj=padj_rest)); wr['sig'] = wr.p_adj < .05
# ---------- runs ----------
bcount = np.array([len(w) for w in run_wells]); ia, ib = np.array(list(itertools.combinations(range(A), 2))).T
def run_stats(wr_):
    m = np.vstack([M[wr_ == i].mean(0) if (wr_ == i).any() else np.zeros(K) for i in range(A)]); b = np.bincount(wr_, minlength=A)
    Tpair = b[ia] * b[ib] / np.maximum(b[ia] + b[ib], 1) * sq(m[ia] - m[ib]) / MSB; tot = M.sum(0); oth = (tot - b[:, None] * m) / np.maximum(B - b, 1)[:, None]
    Tone = b * (B - b) / B * sq(m - oth) / MSB; return m, Tpair, Tone
m_obs, Tpair_o, Tone_o = run_stats(well_run); mx_pair = np.zeros(NP); mx_one = np.zeros(NP)
for b in range(NP): _, Tp, To = run_stats(rng.permutation(well_run)); mx_pair[b] = Tp.max(); mx_one[b] = To.max()
rp = pd.DataFrame(dict(run_a=[runs[i] for i in ia], run_b=[runs[i] for i in ib], wells_a=bcount[ia], wells_b=bcount[ib], W2_distance=np.sqrt(sq(m_obs[ia] - m_obs[ib])), T=Tpair_o,
                       p_adj=[(1 + (mx_pair >= t).sum()) / (1 + NP) for t in Tpair_o])); rp['sig'] = rp.p_adj < .05
oth_o = (M.sum(0) - bcount[:, None] * m_obs) / (B - bcount)[:, None]
rr = pd.DataFrame(dict(run=runs, wells=bcount, W2_to_other_runs=np.sqrt(sq(m_obs - oth_o)), T=Tone_o, p_adj=[(1 + (mx_one >= t).sum()) / (1 + NP) for t in Tone_o])); rr['sig'] = rr.p_adj < .05
for n_, d_ in [('well_pairs', wp), ('well_vs_rest_of_run', wr), ('run_pairs', rp), ('run_vs_rest', rr)]: d_.to_csv(f'{PH}/{n_}.csv', index=False)
pd.set_option('display.width', 200)
print(f'== {tp}: sig well pairs {int(wp.sig.sum())}/{len(wp)}; sig wells vs rest {int(wr.sig.sum())}/{int(wr.p_adj.notna().sum())}; sig run pairs {int(rp.sig.sum())}/{len(rp)}; sig runs vs rest {int(rr.sig.sum())}/{A}')
print(wp[wp.sig].round(4).to_string(index=False)); print(wr[wr.sig].round(4).to_string(index=False)); print(rp[rp.sig].round(4).to_string(index=False)); print(rr.round(4).to_string(index=False))
# ---------- figure ----------
rc = {r: plt.get_cmap('tab10')(i) for i, r in enumerate(runs)}; SIG = '#c0392b'
fig = plt.figure(figsize=(21, 11)); gs = fig.add_gridspec(2, 7, height_ratios=[1.1, 1], hspace=.6, wspace=.55)
ax = fig.add_subplot(gs[0, 0:3]); D = np.zeros((A, A)); P_ = np.ones((A, A))
for k in range(len(ia)): D[ia[k], ib[k]] = D[ib[k], ia[k]] = rp.W2_distance[k]; P_[ia[k], ib[k]] = P_[ib[k], ia[k]] = rp.p_adj[k]
im = ax.imshow(D, cmap='magma_r'); ax.set_xticks(range(A)); ax.set_yticks(range(A)); ax.set_xticklabels([f'{r}\n({b} w)' for r, b in zip(runs, bcount)], fontsize=8); ax.set_yticklabels(runs, fontsize=9)
for i in range(A):
    for j in range(A):
        if i != j: ax.text(j, i, f'{D[i, j]:.2f}' + ('*' if P_[i, j] < .05 else ''), ha='center', va='center', fontsize=8.5, color='white' if D[i, j] > D.max() * .55 else 'black', fontweight='bold' if P_[i, j] < .05 else 'normal')
plt.colorbar(im, ax=ax, fraction=.04, pad=.03); ax.set_title('Run vs run: W2 distance between run barycenters (log10 ISI)\n* = max-T adjusted p < 0.05', fontsize=10)
ax = fig.add_subplot(gs[0, 3:5]); ax.bar(range(A), rr.W2_to_other_runs, color=[SIG if s else '#b9b9b9' for s in rr.sig]); ax.set_xticks(range(A)); ax.set_xticklabels([f'{r}\n({b} w)' for r, b in zip(runs, bcount)], fontsize=8)
for i, (p, s) in enumerate(zip(rr.p_adj, rr.sig)): ax.text(i, rr.W2_to_other_runs[i], f'p={p:.3f}' if p >= .001 else 'p<.001', ha='center', va='bottom', fontsize=7.5, color=SIG if s else '#777')
ax.set_ylabel('W2 distance to the other runs'); ax.set_title('Each run vs the other runs', fontsize=10.5)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
ax = fig.add_subplot(gs[0, 5:7]); ax.axis('off')
txt = f'{tp.upper()}  (DIV {sorted(set(DIVS[tp].values()))}): Frechet post hoc\n\nRUN level (wells are the units, family = all runs):\n'
txt += f'  significant run pairs: ' + (', '.join(f'{a}-{b}' for a, b in zip(rp[rp.sig].run_a, rp[rp.sig].run_b)) or 'none') + '\n  runs differing from the rest: ' + (', '.join(rr[rr.sig].run) or 'none')
txt += '\n\nWELL level (electrodes are the units, family = wells of one run):\n'
for i, r in enumerate(runs):
    s = wp[(wp.run == r) & wp.sig]; v = wr[(wr.run == r) & wr.sig]
    txt += f'  run {r} ({bcount[i]} wells): ' + ('single well, no comparison' if bcount[i] < 2 else ('; '.join(f'w{a}≠w{b}' for a, b in zip(s.well_a, s.well_b)) or 'no pair differs')) + (f'   | outlying vs rest: w' + ','.join(v.well) if len(v) else '') + '\n'
txt += '\nCaveat: runs 22-24 were scanned 30 s, 42/43/45/AR 60 s; full-range ISIs are capped by scan length,\nso run-level differences can be a scan-length artefact. Well-vs-well within a run is not affected.'
ax.text(0, 1, txt, va='top', fontsize=8.4, family='monospace')
for i, r in enumerate(runs):
    axi = fig.add_subplot(gs[1, i]); sub = wr[wr.run == r]; xs = np.arange(len(sub))
    axi.bar(xs, sub.W2_to_rest_of_run.fillna(0), color=[SIG if s else rc[r] for s in sub.sig], alpha=.85); axi.set_xticks(xs); axi.set_xticklabels([f'w{w}' for w in sub.well], fontsize=8)
    for x, (p, s) in enumerate(zip(sub.p_adj, sub.sig)):
        if not np.isnan(p): axi.text(x, sub.W2_to_rest_of_run.iloc[x], '*' if s else '', ha='center', va='bottom', fontsize=14, color=SIG)
    axi.set_title(f'run {r}' + ('  (omnibus p ' + ('<0.001' if R["per"].set_index("run").loc[r, "p"] < .001 else f'= {R["per"].set_index("run").loc[r, "p"]:.3f}') + ')' if r in set(R['per'].run) and not np.isnan(R['per'].set_index('run').loc[r, 'p']) else '  (1 well)'), fontsize=8.5, color=rc[r], fontweight='bold')
    if bcount[i] < 2: axi.text(.5, .5, 'single well:\nno within-run\ncomparison', transform=axi.transAxes, ha='center', va='center', fontsize=9, color='#777'); axi.set_yticks([]); axi.set_xticks([])
    if i == 0: axi.set_ylabel('W2 distance of well to\nrest of its run')
    for s in ('top', 'right'): axi.spines[s].set_visible(False)
fig.suptitle(f'{tp}: nested Frechet ANOVA post hoc — red * = well differs from the rest of its run (max-T adjusted p < 0.05)', y=.98, fontsize=12)
fig.savefig(f'{PH}/posthoc_{tp}.png', dpi=150, bbox_inches='tight', facecolor='white')
