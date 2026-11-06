"""Radial tree: root (grand mean) -> runs (run mean) -> wells (well mean) -> electrodes (bar = mean ISI ms), + variance decomposition + formulas.
Usage: 04_tree_figure.py early|mid|late.  Node numbers are means of ln(mean ISI [ms]) (the ANOVA scale); the ms equivalent (geometric mean) is printed under well nodes."""
import pickle
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from common import *
tp, OUT = setup()
R = pickle.load(open(f'{OUT}/_res.pkl', 'rb')); th, res, vc, pt = R['th'].copy(), R['res'], R['vc'], R['pt']; D_MIN, MIN_M = R['D_MIN'], R['MIN_M']
th['y'] = np.log(th.mean_isi_ms)
INK, MUT = '#1a1a1a', '#8a8a8a'
LV = {'run': '#2f6db3', 'well': '#e08a1e', 'electrode': '#b9b9b9'}
runs = sorted(th.run.unique(), key=lambda r: (not r.isdigit(), int(r) if r.isdigit() else 0))
rc = {r: plt.get_cmap('tab10')(i) for i, r in enumerate(runs)}
grand = th.y.mean(); run_mean = th.groupby('run').y.mean(); well_mean = th.groupby(['run', 'well']).y.mean()
r_run, r_well, r_el = 1.1, 2.5, 3.4; GAP = 0.06; GAPW = 0.008; R_LAB = 4.3
n = len(th); LABELS = n <= 600                                    # ring of per-electrode ms numbers only legible when few electrodes
msz = th.groupby(['run', 'well']).size(); weight = 0.25 * msz / msz.mean() + 0.75          # wedge per well: mostly equal share so nodes never overlap
unit = (2 * np.pi - GAP * len(runs) - GAPW * (len(msz) - len(runs))) / weight.sum(); bar_scale = 0.8 / th.mean_isi_ms.max(); blw = 1.1 if n < 700 else 0.7 if n < 1800 else 0.45
pol = lambda t, r: (r * np.cos(t), r * np.sin(t))
divs = sorted(set(DIVS[tp].values()))

fig = plt.figure(figsize=(27, 15.5)); gs = GridSpec(2, 2, width_ratios=[1.75, 1], height_ratios=[0.8, 1.2], wspace=0.04, hspace=0.12)
ax = fig.add_subplot(gs[:, 0]); ax.set_aspect('equal'); ax.axis('off')
t = 0.
for r in runs:
    t0 = t; wa_list = []
    for w in sorted(th[th.run == r].well.unique()):
        g = th[(th.run == r) & (th.well == w)].sort_values(['px', 'py']); wd = weight[(r, w)] * unit; ts = t + wd * (np.arange(len(g)) + .5) / len(g); wa = t + wd / 2; t += wd + GAPW
        wa_list.append(wa); x0, y0 = pol(wa, r_well)
        for tt, v in zip(ts, g.mean_isi_ms):
            xa, ya = pol(tt, r_el); ax.plot([x0, xa], [y0, ya], color=LV['electrode'], lw=.2 if n > 600 else .35, alpha=.3 if n > 600 else .6)
            xb, yb = pol(tt, r_el + v * bar_scale); ax.plot([xa, xb], [ya, yb], color=rc[r], lw=blw, solid_capstyle='butt')
            if LABELS:
                deg = np.degrees(tt) % 360; flip = 90 < deg < 270; xl, yl = pol(tt, R_LAB)
                ax.text(xl, yl, f'{v:.0f}', fontsize=5.6, color=rc[r], rotation=deg + (180 if flip else 0), rotation_mode='anchor',
                        ha='right' if flip else 'left', va='center')
        ax.plot(*pol(wa, r_well), 'o', ms=23, color='white', mec=rc[r], mew=1.6, zorder=5)
        wm = well_mean[(r, w)]
        ax.text(*pol(wa, r_well), f'{wm:.2f}\n{np.exp(wm):.0f} ms', fontsize=5.2, ha='center', va='center', zorder=6, color=INK, linespacing=1.15)
        ax.text(*pol(wa, r_well + .27), f'w{w}', fontsize=6.5, ha='center', va='center', color=MUT)
    ra = (t0 + t) / 2; xr, yr = pol(ra, r_run)
    for wa in wa_list: ax.plot([xr, pol(wa, r_well)[0]], [yr, pol(wa, r_well)[1]], color=LV['well'], lw=1.6, alpha=.9)
    ax.plot([0, xr], [0, yr], color=LV['run'], lw=2.2, alpha=.9)
    ax.plot(xr, yr, 's', ms=36, color='white', mec=rc[r], mew=2.2, zorder=6)
    ax.text(xr, yr, f'{run_mean[r]:.2f}\n{np.exp(run_mean[r]):.0f} ms', fontsize=8, fontweight='bold', ha='center', va='center', zorder=7, color=INK, linespacing=1.15)
    ax.text(*pol(ra, R_LAB + .75), f'run {r}', fontsize=11, fontweight='bold', color=rc[r], ha='center', va='center')
    t += GAP - GAPW
ax.plot(0, 0, 'o', ms=52, color=INK, zorder=7)
ax.text(0, 0.17, 'grand mean', color='white', fontsize=6.5, ha='center', va='center', zorder=8)
ax.text(0, -0.02, f'{grand:.2f}', color='white', fontsize=10, fontweight='bold', ha='center', va='center', zorder=8)
ax.text(0, -0.19, f'{np.exp(grand):.0f} ms', color='white', fontsize=7, ha='center', va='center', zorder=8)
for rv in (50, 100, 150):
    ax.add_patch(plt.Circle((0, 0), r_el + rv * bar_scale, fill=False, ls=':', lw=.5, color=MUT)); ax.text(0.05, r_el + rv * bar_scale, f'{rv} ms', fontsize=6.5, color=MUT, va='bottom')
ax.set_xlim(-5.7, 5.7); ax.set_ylim(-5.7, 5.7)
ax.legend(handles=[Line2D([], [], color=LV['run'], lw=3, label='run level:  root → run  (run mean vs grand mean)'),
                   Line2D([], [], color=LV['well'], lw=3, label='well level:  run → well  (well mean vs run mean)'),
                   Line2D([], [], color=LV['electrode'], lw=3, label='electrode level:  well → electrode  (electrode vs well mean)')],
          loc='lower left', fontsize=9, frameon=False, bbox_to_anchor=(0, -.01))
ax.set_title(f'{tp.upper()} time point (DIV {divs}), vehicle control.  Node numbers = mean of ln(mean ISI [ms]) over all electrodes below the node; '
             f'second line = exp(mean) in ms (geometric mean).\nOuter bars' + (' and ring of numbers' if LABELS else '') + f' = each electrode mean ISI in ms (<200 ms)'
             + ('' if LABELS else f'  (per-electrode numbers omitted: {n} electrodes, in electrodes_thinned.csv)')
             + f'.\nElectrodes thinned to ≥{D_MIN:.0f} µm apart (well_counts.csv: electrodes kept per well).', fontsize=10, loc='left')

# ---- variance bars ----
ax2 = fig.add_subplot(gs[0, 1]); order = ['ANOVA thinned', 'REML thinned', 'ANOVA all electrodes', 'REML all electrodes']; left = np.zeros(4)
for k_ in ['run', 'well', 'electrode']:
    vals = np.array([vc.set_index('method').loc[o, f'prop_{k_}'] * 100 for o in order])
    ax2.barh(range(4), vals, left=left, color=LV[k_], label={'run': 'run  σ²α', 'well': 'well(run)  σ²β', 'electrode': 'electrode  σ²'}[k_], height=.6)
    for i, (l, v) in enumerate(zip(left, vals)):
        if v > 4: ax2.text(l + v / 2, i, f'{v:.1f}%', ha='center', va='center', fontsize=8.5)
        else: ax2.text(l + v + .6, i - .38 if k_ == 'run' else i + .38, f'{v:.1f}%', fontsize=7.5, color=LV[k_], va='center')
    left += vals
ax2.set_yticks(range(4)); ax2.set_yticklabels(order, fontsize=9); ax2.invert_yaxis(); ax2.set_xlim(0, 100)
ax2.set_xlabel('% of total variance of ln(mean ISI)'); ax2.legend(ncol=3, fontsize=8, frameon=False, loc='lower center', bbox_to_anchor=(.5, 1.0))
for s in ('top', 'right'): ax2.spines[s].set_visible(False)

# ---- how variance is computed, typeset ----
ax3 = fig.add_subplot(gs[1, 1]); ax3.axis('off'); ax3.set_xlim(0, 1); ax3.set_ylim(0, 1)
T = res['tab'].set_index('source'); k1, k2, k3 = res['k']; c = res['comp']; tot = sum(c.values()); N = res['N']
tr = lambda k: f"{c[k]/tot*100:.1f}%"
fp = lambda p: r'<0.0001' if p < 1e-4 else f'{p:.4f}' if p < .01 else f'{p:.3f}'
pe = lambda p: ('' if p < 1e-4 else '=') + fp(p)
def block(y, color, title, rows):
    ax3.add_patch(plt.Rectangle((0, y - .215), .012, .245, color=color, clip_on=False))
    ax3.text(.03, y + .02, title, fontsize=10.5, fontweight='bold', va='center', color=INK)
    for i, r_ in enumerate(rows): ax3.text(.04, y - .045 - i * .075, r_, fontsize=11.5, va='center', color=INK)
block(.93, LV['run'], 'RUN (blue edges): how far run means sit from the grand mean', [
    rf"$SS_{{run}}=\sum_{{i}} n_{{i}}\,(\bar{{y}}_{{i\cdot\cdot}}-\bar{{y}}_{{\cdot\cdot\cdot}})^2={T.loc['Run','SS']:.2f}$,   $df=a-1={int(T.loc['Run','df'])}$,   $MS_{{run}}={T.loc['Run','MS']:.3f}$",
    rf"$\hat\sigma^2_{{run}}=\dfrac{{MS_{{run}}-MS_{{e}}-k_2\,\hat\sigma^2_{{well}}}}{{k_3}}={c['run']:.4f}$  →  {tr('run')} of total",
    rf"$F=\dfrac{{MS_{{run}}}}{{MS_{{well(run)}}}}={T.loc['Run','F']:.2f}$:   F-test $p{pe(T.loc['Run','p'])}$,   permutation $p{pe(pt['pA'])}$"])
block(.66, LV['well'], 'WELL within run (orange edges): how far well means sit from their run mean', [
    rf"$SS_{{well}}=\sum_{{i}}\sum_{{j}} n_{{ij}}\,(\bar{{y}}_{{ij\cdot}}-\bar{{y}}_{{i\cdot\cdot}})^2={T.loc['Well(Run)','SS']:.2f}$,   $df={int(T.loc['Well(Run)','df'])}$,   $MS_{{well}}={T.loc['Well(Run)','MS']:.3f}$",
    rf"$\hat\sigma^2_{{well}}=\dfrac{{MS_{{well}}-MS_{{e}}}}{{k_1}}={c['well']:.4f}$  →  {tr('well')} of total",
    rf"$F=\dfrac{{MS_{{well}}}}{{MS_{{e}}}}={T.loc['Well(Run)','F']:.2f}$:   F-test $p{pe(T.loc['Well(Run)','p'])}$,   permutation $p{pe(pt['pB'])}$"])
block(.39, LV['electrode'], 'ELECTRODE (grey edges): how far electrodes sit from their well mean', [
    rf"$SS_{{e}}=\sum_{{i}}\sum_{{j}}\sum_{{k}} (y_{{ijk}}-\bar{{y}}_{{ij\cdot}})^2={T.loc['Electrode (resid)','SS']:.1f}$,   $df=N-B={int(T.loc['Electrode (resid)','df'])}$",
    rf"$\hat\sigma^2_{{e}}=MS_{{e}}={c['electrode']:.4f}$  →  {tr('electrode')} of total",
    r"not testable: one value per electrode, so it is the residual"])
ax3.text(0, .075, rf"$\hat\sigma^2_{{total}}=\hat\sigma^2_{{run}}+\hat\sigma^2_{{well}}+\hat\sigma^2_{{e}}={tot:.4f}$   ($N={N}$ electrodes, $B={res['B']}$ wells, $a={res['a']}$ runs, unbalanced)",
         fontsize=10.5, va='center', color=INK)
ax3.text(0, -.025, rf"$\bar{{y}}_{{ij\cdot}}$ = well mean, $\bar{{y}}_{{i\cdot\cdot}}$ = run mean (over all its electrodes), $\bar{{y}}_{{\cdot\cdot\cdot}}$ = grand mean.""" + "\n" +
         rf"Searle coefficients $k_1={k1:.2f},\ k_2={k2:.2f},\ k_3={k3:.2f}$ (equal $n_{{ij}}$: $n,\ n,\ bn$).  F-test = parametric (assumes independent normal errors); permutation = distribution-free.  Bars at top = $\hat\sigma^2$ / total.",
         fontsize=8.6, va='center', color=MUT)
fig.savefig(f'{OUT}/nested_anova_tree_{tp}_vc.png', dpi=170, bbox_inches='tight', facecolor='white'); print('saved', tp, n, 'electrodes')
