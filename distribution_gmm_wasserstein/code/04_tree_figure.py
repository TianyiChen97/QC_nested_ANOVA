"""Radial tree of distributions.  root (grand barycenter) -> runs (run barycenters) -> wells (well barycenters) -> electrodes (fitted 2-GMM ISI density).
Usage: 04_tree_figure.py early|mid|late [N_MIN=30] [gmm|emp]
Outer ring: one column per electrode (angle), radius = log10 ISI (ms), colour = that electrode's fitted ISI density (electrodes inside a well sorted by their median ISI, display only).
Node insets: density of the barycenter = distribution whose quantile function is the mean of the quantile functions below the node (KDE of Qbar(u_k), h = 0.08)."""
import pickle
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from scipy.stats import norm
from common import *
tp, OUT0 = setup(); N_MIN = int(sys.argv[2]) if len(sys.argv) > 2 else 30; MODE = sys.argv[3] if len(sys.argv) > 3 else 'gmm'
OUT = f'{OUT0}/nmin{N_MIN}' + ('' if MODE == 'gmm' else '/sensitivity_empirical_quantiles')
R = pickle.load(open(f'{OUT}/_res.pkl', 'rb')); th, Q, runs, wells, o, res_tab, per, comp, pA, pB = R['th'], R['Q'], R['runs'], R['wells'], R['o'], R['tab'], R['per'], R['comp'], R['pA'], R['pB']
th['run'] = th.run.astype(str); th['wid'] = th.run + '_w' + th.well.astype(str); widx = R['widx']; well_run = R['well_run']; k1, k2, k3 = R['k']; vc = R['vc']
INK, MUT = '#1a1a1a', '#8a8a8a'; LV = {'run': '#2f6db3', 'well': '#e08a1e', 'electrode': '#b9b9b9'}
rc = {r: plt.get_cmap('tab10')(i) for i, r in enumerate(runs)}
XL = (0.3, 4.6); XC = np.linspace(*XL, 110); H = 0.08                    # log10 ISI(ms) display range, KDE bandwidth
def gmm_pdf(r): return r.w1 * norm.pdf(XC, r.mu1, r.s1) + r.w2 * norm.pdf(XC, r.mu2, r.s2)
def kde(qv): return norm.pdf(XC[:, None], qv[None, :], H).mean(1)
if MODE == 'gmm': Dens = np.vstack([gmm_pdf(r) for r in th.itertuples()])
else: Dens = np.vstack([kde(q) for q in Q])
med = Q[:, Q.shape[1] // 2]
dW = np.vstack([kde(q) for q in o['M']]); dR = np.vstack([kde(q) for q in o['Mi']]); dG = kde(o['G'])
vmax = np.percentile(Dens, 99.5)
n = len(th); msz = pd.Series(widx).value_counts().sort_index().values
weight = 0.25 * msz / msz.mean() + 0.75; GAP, GAPW = 0.07, 0.008
unit = (2 * np.pi - GAP * len(runs) - GAPW * (len(wells) - len(runs))) / weight.sum()
r_run, r_well, r_in, r_out = 1.35, 3.1, 4.2, 5.3; wrad = {}; pol = lambda t, r: (r * np.cos(t), r * np.sin(t))
divs = sorted(set(DIVS[tp].values()))
fig = plt.figure(figsize=(27, 15.5)); gs = GridSpec(3, 2, width_ratios=[1.75, 1], height_ratios=[0.62, 0.55, 1.2], wspace=0.04, hspace=0.3)
ax = fig.add_subplot(gs[:, 0]); ax.set_aspect('equal'); ax.axis('off'); ax.set_xlim(-6.6, 6.6); ax.set_ylim(-6.6, 6.6)
xe = np.linspace(XL[0] - (XC[1] - XC[0]) / 2, XL[1] + (XC[1] - XC[0]) / 2, len(XC) + 1); re = r_in + (xe - XL[0]) / (XL[1] - XL[0]) * (r_out - r_in)
def inset(x, y, w, h, d, color, ymax, lab=None):
    a = ax.inset_axes([x - w / 2, y - h / 2, w, h], transform=ax.transData, zorder=6); a.set_facecolor((1, 1, 1, .92)); a.fill_between(XC, d, color=color, alpha=.55, lw=0); a.plot(XC, d, color=color, lw=1)
    a.set_xlim(*XL); a.set_ylim(0, ymax); a.set_xticks([]); a.set_yticks([])
    if lab: a.text(.04, .86, lab, transform=a.transAxes, fontsize=6.5, color=color, va='center', fontweight='bold')
    for s in a.spines.values(): s.set_edgecolor(color); s.set_linewidth(1)
t = 0.; run_ang = {}; last_end = 0
for i, r in enumerate(runs):
    t0 = t; wa_list = []
    for j in np.where(well_run == i)[0]:
        wd = weight[j] * unit; idx = np.where(widx == j)[0]; idx = idx[np.argsort(med[idx])]; m = len(idx); te = t + wd * np.arange(m + 1) / m
        TH, RR = np.meshgrid(te, re); ax.pcolormesh(RR * np.cos(TH), RR * np.sin(TH), Dens[idx].T, cmap='magma_r', vmin=0, vmax=vmax, rasterized=True, shading='flat', zorder=2)
        ax.plot(*zip(*[pol(a_, r_out + .06) for a_ in np.linspace(t, t + wd, 40)]), color=rc[r], lw=3.2, solid_capstyle='butt', zorder=3)
        wa = t + wd / 2; wa_list.append(wa); wrad[j] = r_well + (-.22 if len(wa_list) % 2 else .24); x0, y0 = pol(wa, wrad[j])
        for a_ in (t, t + wd): ax.plot([x0, pol(a_, r_in)[0]], [y0, pol(a_, r_in)[1]], color=LV['electrode'], lw=.5, alpha=.7, zorder=1)
        t += wd + GAPW
    ra = (t0 + t) / 2; run_ang[r] = ra; xr, yr = pol(ra, r_run)
    for wa, j in zip(wa_list, np.where(well_run == i)[0]): ax.plot([xr, pol(wa, wrad[j])[0]], [yr, pol(wa, wrad[j])[1]], color=LV['well'], lw=1.6, alpha=.9, zorder=1)
    ax.plot([0, xr], [0, yr], color=LV['run'], lw=2.4, alpha=.9, zorder=1)
    ax.text(*pol(ra, r_out + .6), f'run {r}', fontsize=11.5, fontweight='bold', color=rc[r], ha='center', va='center')
    t += GAP - GAPW; last_end = t
ymW = np.percentile(dW, 99); ymR = dR.max() * 1.05
# node insets (positions recomputed from the same layout)
t = 0.
for i, r in enumerate(runs):
    for j in np.where(well_run == i)[0]:
        wd = weight[j] * unit; wa = t + wd / 2; x0, y0 = pol(wa, wrad[j]); inset(x0, y0, .6, .36, dW[j], rc[r], ymW, 'w' + wells[j].split('_w')[1]); t += wd + GAPW
    xr, yr = pol(run_ang[r], r_run); inset(xr, yr, 1.0, .66, dR[i], rc[r], ymR); t += GAP - GAPW
inset(0, 0, 1.5, .95, dG, INK, ymR)
for xv, lab in [(1, '10 ms'), (2, '100 ms'), (3, '1 s'), (4, '10 s')]:  # radial axis of the ring, drawn in the last gap
    rr = r_in + (xv - XL[0]) / (XL[1] - XL[0]) * (r_out - r_in); ang = last_end - GAP / 2 + GAPW / 2
    ax.plot(*zip(*[pol(a_, rr) for a_ in np.linspace(ang - .035, ang + .035, 5)]), color=INK, lw=.9, zorder=4); ax.text(pol(ang, rr)[0], pol(ang, rr)[1] + .1, lab, fontsize=6.5, color=INK, ha='center', va='bottom', rotation=90)
ax.legend(handles=[Line2D([], [], color=LV['run'], lw=3, label='run level:  root → run  (run barycenter vs grand barycenter)'),
                   Line2D([], [], color=LV['well'], lw=3, label='well level:  run → well  (well barycenter vs run barycenter)'),
                   Line2D([], [], color=LV['electrode'], lw=3, label='electrode level:  well → electrode  (electrode distribution vs well barycenter)')], loc='lower left', fontsize=9, frameon=False)
sm = plt.cm.ScalarMappable(cmap='magma_r', norm=plt.Normalize(0, vmax)); cax = ax.inset_axes([.80, .035, .17, .012]); cb = fig.colorbar(sm, cax=cax, orientation='horizontal'); cb.set_label('electrode ISI density (per log10 unit)', fontsize=8); cb.ax.tick_params(labelsize=7)
ax.set_title(f'{tp.upper()} time point (DIV {divs}), vehicle control.  Each electrode = its ISI distribution: ' + ('2-component GMM fitted to log10(ISI), shown as a density column in the outer ring' if MODE == 'gmm' else 'empirical (KDE of its ISIs)') +
             f'; radius = ISI (log scale), electrodes sorted by median within each well.\nNode insets = barycenter density (distribution with the mean quantile function of everything below the node).  '
             f'\nElectrodes need ≥{N_MIN} ISIs; thinned to ≥{D_BY_TP[tp]:.0f} µm apart.  N = {n} electrodes, {len(wells)} wells, {len(runs)} runs.', fontsize=10, loc='left')

# ---- variance bars ----
ax2 = fig.add_subplot(gs[0, 1]); rows = vc.method.tolist(); lab = ['Fréchet ANOVA, thinned (tested)', 'all eligible electrodes, empirical Q (descriptive)']; left = np.zeros(2)
for k_ in ['run', 'well', 'electrode']:
    vals = np.array([vc.iloc[i][f'prop_{k_}'] * 100 for i in range(2)])
    ax2.barh(range(2), vals, left=left, color=LV[k_], height=.55, label={'run': 'run  σ²α', 'well': 'well(run)  σ²β', 'electrode': 'electrode  σ²'}[k_])
    for i, (l, v) in enumerate(zip(left, vals)):
        if v > 4: ax2.text(l + v / 2, i, f'{v:.1f}%', ha='center', va='center', fontsize=9)
        else: ax2.text(l + v + .6, i - .36 if k_ == 'run' else i + .36, f'{v:.1f}%', fontsize=7.5, color=LV[k_], va='center')
    left += vals
ax2.set_yticks(range(2)); ax2.set_yticklabels(lab, fontsize=9); ax2.invert_yaxis(); ax2.set_xlim(0, 100); ax2.set_xlabel('% of total Fréchet variance of the ISI distribution (W2², log10 ms)')
ax2.legend(ncol=3, fontsize=8, frameon=False, loc='lower center', bbox_to_anchor=(.5, 1.0))
for s in ('top', 'right'): ax2.spines[s].set_visible(False)
# ---- variance profile over quantile level ----
ax4 = fig.add_subplot(gs[1, 1]); u = U
for k_, v in [('run', R['ca']), ('well', R['cb']), ('electrode', R['ce'])]: ax4.plot(u, v, color=LV[k_] if k_ != 'electrode' else '#7a7a7a', lw=2, label={'run': 'run', 'well': 'well(run)', 'electrode': 'electrode'}[k_])
ax4.set_yscale('log'); ax4.set_xlabel('quantile level u of the ISI distribution (0 = shortest ISIs, 1 = longest)'); ax4.set_ylabel('variance at level u\n(log10 ISI)²'); ax4.legend(fontsize=8, frameon=False, ncol=3)
for s in ('top', 'right'): ax4.spines[s].set_visible(False)
# ---- formulas ----
ax3 = fig.add_subplot(gs[2, 1]); ax3.axis('off'); ax3.set_xlim(0, 1); ax3.set_ylim(0, 1)
T = res_tab.set_index('source'); c = comp; tot = sum(c.values()); N = len(th); tr = lambda k: f"{c[k]/tot*100:.1f}%"
fp = lambda p: r'<0.0002' if p <= 1.0001 / 5001 else f'{p:.4f}'; pe = lambda p: ('' if p <= 1.0001 / 5001 else '=') + fp(p)
def block(y, color, title, rows_):
    ax3.add_patch(plt.Rectangle((0, y - .215), .012, .245, color=color, clip_on=False)); ax3.text(.03, y + .02, title, fontsize=10.5, fontweight='bold', va='center', color=INK)
    for i, r_ in enumerate(rows_): ax3.text(.04, y - .045 - i * .075, r_, fontsize=11.2, va='center', color=INK)
block(.95, LV['run'], 'RUN (blue edges): how far run barycenters sit from the grand barycenter', [
    rf"$SS_{{run}}=\sum_{{i}} n_{{i}}\,\Vert \bar{{Q}}_{{i\cdot\cdot}}-\bar{{Q}}_{{\cdot\cdot\cdot}}\Vert ^2_{{L^2}}={T.loc['Run','SS']:.2f}$,   $df={int(T.loc['Run','df'])}$,   $MS={T.loc['Run','MS']:.3f}$",
    rf"$\hat\sigma^2_{{run}}=\dfrac{{MS_{{run}}-MS_{{e}}-k_2\,\hat\sigma^2_{{well}}}}{{k_3}}={c['run']:.4f}$  →  {tr('run')} of total",
    rf"$F^{{*}}=\dfrac{{MS_{{run}}}}{{MS_{{well(run)}}}}={T.loc['Run','F_stat']:.2f}$,   permutation $p{pe(pA)}$"])
block(.68, LV['well'], 'WELL within run (orange edges): how far well barycenters sit from their run barycenter', [
    rf"$SS_{{well}}=\sum_{{i,j}} n_{{ij}}\,\Vert \bar{{Q}}_{{ij\cdot}}-\bar{{Q}}_{{i\cdot\cdot}}\Vert ^2_{{L^2}}={T.loc['Well(Run)','SS']:.2f}$,   $df={int(T.loc['Well(Run)','df'])}$,   $MS={T.loc['Well(Run)','MS']:.3f}$",
    rf"$\hat\sigma^2_{{well}}=\dfrac{{MS_{{well}}-MS_{{e}}}}{{k_1}}={c['well']:.4f}$  →  {tr('well')} of total",
    rf"$F^{{*}}=\dfrac{{MS_{{well}}}}{{MS_{{e}}}}={T.loc['Well(Run)','F_stat']:.2f}$,   permutation $p{pe(pB)}$"])
block(.41, LV['electrode'], 'ELECTRODE (grey edges): how far each electrode distribution sits from its well barycenter', [
    rf"$SS_{{e}}=\sum_{{i,j,k}} \Vert Q_{{ijk}}-\bar{{Q}}_{{ij\cdot}}\Vert ^2_{{L^2}}={T.loc['Electrode (resid)','SS']:.2f}$,   $df=N-B={int(T.loc['Electrode (resid)','df'])}$",
    rf"$\hat\sigma^2_{{e}}=MS_{{e}}={c['electrode']:.4f}$  →  {tr('electrode')} of total", r"$W_2^2$ distance between ISI distributions = $\Vert Q_1-Q_2\Vert ^2_{L^2}$ (quantile functions of log10 ISI)"])
ax3.text(0, .105, rf"$\hat\sigma^2_{{total}}=\hat\sigma^2_{{run}}+\hat\sigma^2_{{well}}+\hat\sigma^2_{{e}}={tot:.4f}$   ($N={N}$, $B={len(wells)}$, $a={len(runs)}$, unbalanced);   $\Vert f\Vert ^2_{{L^2}}=\int_0^1 f(u)^2du$",
         fontsize=10.2, va='center', color=INK)
ax3.text(0, .02, rf"$\bar{{Q}}_{{ij\cdot}}$ = well barycenter = mean quantile function of its electrodes; $\bar{{Q}}_{{i\cdot\cdot}}$ run; $\bar{{Q}}_{{\cdot\cdot\cdot}}$ grand.  Searle $k_1={k1:.2f},\ k_2={k2:.2f},\ k_3={k3:.2f}$." + "\n" +
         r"$F^{*}$ has no F reference for functional data: p-values are permutation only (run: permute run labels of whole wells; well: permute electrodes among wells within run).",
         fontsize=8.3, va='center', color=MUT)
fig.savefig(f'{OUT}/frechet_tree_{tp}_vc.png', dpi=150, bbox_inches='tight', facecolor='white'); print('saved', tp, n)
