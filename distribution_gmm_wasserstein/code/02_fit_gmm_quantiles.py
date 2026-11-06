"""Per-electrode 2-component Gaussian mixture on log10(ISI ms) -> quantile function on the grid U.   Usage: 02_fit_gmm_quantiles.py early|mid|late [N_MIN=30] [D_MIN_um]
1. eligible electrodes = active electrodes with >= N_MIN ISIs.   2. thin each well to eligible electrodes >= D_MIN um apart (positions only; random order, seed 0).
3. keep wells with >= MIN_M thinned electrodes.   4. fit GMM (K=2, full var, reg_covar=1e-3 -> sd floor 0.03 log10 units, 3 inits), components ordered by mean.
5. mixture CDF -> quantile function Q(u_k), k=1..K, u_k=(k-.5)/K  (label-switch invariant).   Also saves empirical quantiles (no model) for sensitivity."""
import os, sys, pickle
import numpy as np, pandas as pd
from scipy.stats import norm
from sklearn.mixture import GaussianMixture
from joblib import Parallel, delayed
from common import *
tp, OUT0 = setup()
N_MIN = int(sys.argv[2]) if len(sys.argv) > 2 else 30
D_MIN = float(sys.argv[3]) if len(sys.argv) > 3 else D_BY_TP[tp]
MIN_M = 10; rng = np.random.default_rng(0)
OUT = f'{OUT0}/nmin{N_MIN}'; os.makedirs(OUT, exist_ok=True)
meta = pd.read_csv(f'{OUT0}/electrodes_isi_meta.csv'); isis = pickle.load(open(f'{OUT0}/isis.pkl', 'rb'))
meta['run'] = meta.run.astype(str); meta['wid'] = meta.run + '_w' + meta.well.astype(str)
elig = meta[meta.n_isi >= N_MIN].reset_index(drop=True)
lx = lambda key: np.log10(np.clip(isis[key].astype(float), 0.1, None))

def thin(g):
    xy = g[['px', 'py']].values; kept = []
    for i in rng.permutation(len(g)):
        if all(np.hypot(*(xy[i] - xy[k])) >= D_MIN for k in kept): kept.append(i)
    return g.iloc[sorted(kept)]
thinned = pd.concat([thin(g) for _, g in elig.groupby('wid')]) if len(elig) else elig
cnt = pd.DataFrame({'eligible': elig.groupby('wid').size(), 'm_thinned': thinned.groupby('wid').size()}).fillna(0).astype(int)
# all 41-42 wells listed, including those with no eligible electrodes
allw = sorted({f'{r}_w{w}' for r in DIVS[tp] for w in range(1, 7)}); cnt = cnt.reindex(allw).fillna(0).astype(int); cnt.to_csv(f'{OUT}/well_counts.csv')
keep = cnt.index[cnt.m_thinned >= MIN_M]; th = thinned[thinned.wid.isin(keep)].reset_index(drop=True)
print(f'{tp}: N_MIN={N_MIN} d={D_MIN:.0f}um | eligible electrodes {len(elig)} | wells kept (>= {MIN_M} thinned) {len(keep)}/{len(cnt)} | thinned electrodes kept {len(th)}')
if len(keep) < 3 or th.run.nunique() < 2: print('NOT ENOUGH DATA for the nested analysis'); open(f'{OUT}/NOT_ENOUGH_DATA.txt', 'w').write(f'{len(keep)} wells, {th.run.nunique() if len(th) else 0} runs'); sys.exit(0)

def fit(key):
    x = lx(key)[:, None]; n = len(x)
    g2 = GaussianMixture(2, n_init=3, reg_covar=1e-3, random_state=0).fit(x)
    o = np.argsort(g2.means_[:, 0]); w = g2.weights_[o]; mu = g2.means_[o, 0]; sd = np.sqrt(g2.covariances_[o, 0, 0])
    return dict(key=key, w1=w[0], mu1=mu[0], s1=sd[0], w2=w[1], mu2=mu[1], s2=sd[1], loglik=g2.score(x) * n, bic2=g2.bic(x),
                bic1=GaussianMixture(1, reg_covar=1e-3).fit(x).bic(x))
fits = pd.DataFrame(Parallel(n_jobs=-1)(delayed(fit)(k) for k in th.key)); th = th.merge(fits, on='key')

xs = np.linspace(-1.5, 5.0, 6501)
def mix_quantiles(r):
    F = r.w1 * norm.cdf(xs, r.mu1, r.s1) + r.w2 * norm.cdf(xs, r.mu2, r.s2); return np.interp(U, F, xs)
Qg = np.vstack([mix_quantiles(r) for r in th.itertuples()])
Qe = np.vstack([np.quantile(lx(k), U) for k in th.key])                       # empirical quantile functions (sensitivity, no model)
th['w2_fit_vs_emp'] = np.sqrt(((Qg - Qe) ** 2).mean(1))                       # 2-Wasserstein distance fit vs data (log10 units)
# BIC check on a random subset of 300 electrodes: K=2 vs K=3
sub = th.sample(min(300, len(th)), random_state=0)
def bic3(key): x = lx(key)[:, None]; return GaussianMixture(3, n_init=3, reg_covar=1e-3, random_state=0).fit(x).bic(x)
b3 = pd.Series(Parallel(n_jobs=-1)(delayed(bic3)(k) for k in sub.key), index=sub.index)
th['bic3'] = np.nan; th.loc[sub.index, 'bic3'] = b3
th.to_csv(f'{OUT}/electrode_fits.csv', index=False); np.save(f'{OUT}/Q_gmm.npy', Qg); np.save(f'{OUT}/Q_emp.npy', Qe)
# empirical quantiles of ALL eligible electrodes (descriptive variance components, no thinning)
np.save(f'{OUT}/Q_emp_all.npy', np.vstack([np.quantile(lx(k), U) for k in elig.key])); elig.to_csv(f'{OUT}/eligible_all_meta.csv', index=False)
d2 = th.bic2 < th.bic1; s = th.dropna(subset=['bic3'])
print(f'  GMM K=2 preferred over K=1 by BIC: {d2.mean()*100:.1f}% of thinned electrodes; vs K=3 (n={len(s)}): K=2 better in {(s.bic2 < s.bic3).mean()*100:.1f}%')
print(f'  W2(fit, empirical) median {th.w2_fit_vs_emp.median():.3f} log10 units; median n_isi {int(th.n_isi.median())}')
