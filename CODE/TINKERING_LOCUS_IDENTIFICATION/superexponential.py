#!/usr/bin/env python3
"""
Fit nearest-neighbor distance data to:
  1. Single exponential:  f(x) = lam * exp(-lam * x)
  2. Mixture of two exponentials: f(x) = pi1*lam1*exp(-lam1*x) + pi2*lam2*exp(-lam2*x)

Fits on raw distances (not histogram), plots both fits atop a histogram.
Reports bootstrapped KS p-values for both models, AIC, BIC, and a
parametric-bootstrap likelihood ratio test.

Usage:
    python fit_distances.py distances.txt
    python fit_distances.py distances.txt --bins 100 --output fit.png
    python fit_distances.py distances.txt --n-boot 9999 --seed 42 --plot-percentile 99.5
"""

__version__ = '2.1.0'

import sys
import json
import argparse
import multiprocessing
import numpy as np
import scipy
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import kstest


# ---------------------------------------------------------------------------
# Single exponential
# ---------------------------------------------------------------------------

def single_nll(params, distances):
    log_lam, = params
    lam = np.exp(log_lam)
    return -np.sum(np.log(lam) - lam * distances)


def fit_single(distances):
    lam0 = 1.0 / np.mean(distances)
    res = minimize(
        single_nll,
        x0=[np.log(lam0)],
        args=(distances,),
        method='Nelder-Mead',
        options={'maxiter': 50000, 'xatol': 1e-10, 'fatol': 1e-10}
    )
    return res


def single_cdf(x, lam):
    return 1.0 - np.exp(-lam * x)


def sample_single(lam, n, rng):
    """Draw n samples from Exp(lam)."""
    return rng.exponential(scale=1.0 / lam, size=n)


# ---------------------------------------------------------------------------
# Mixture of two exponentials
# ---------------------------------------------------------------------------

def mixture_nll(params, distances):
    pi1, log_lam1, log_lam2 = params
    if not (0 < pi1 < 1):
        return 1e12
    lam1 = np.exp(log_lam1)
    lam2 = np.exp(log_lam2)
    pi2  = 1.0 - pi1
    pdf  = pi1 * lam1 * np.exp(-lam1 * distances) + pi2 * lam2 * np.exp(-lam2 * distances)
    # Guard against log(0)
    pdf  = np.maximum(pdf, 1e-300)
    return -np.sum(np.log(pdf))


def fit_mixture(distances):
    lam0 = 1.0 / np.mean(distances)
    lam_hi_candidates = [lam0 * 5,   lam0 * 20,  lam0 * 50]
    lam_lo_candidates = [lam0 * 0.1, lam0 * 0.3, lam0 * 0.5]
    pi1_candidates    = [0.05, 0.1, 0.2]

    best_result = None
    best_nll    = np.inf

    for pi1_0 in pi1_candidates:
        for lam_hi in lam_hi_candidates:
            for lam_lo in lam_lo_candidates:
                if lam_hi <= lam_lo:
                    continue
                x0 = [pi1_0, np.log(lam_hi), np.log(lam_lo)]
                res = minimize(
                    mixture_nll,
                    x0=x0,
                    args=(distances,),
                    method='Nelder-Mead',
                    options={'maxiter': 50000, 'xatol': 1e-8, 'fatol': 1e-8}
                )
                if res.fun < best_nll:
                    best_nll    = res.fun
                    best_result = res

    return best_result


def fit_mixture_fast(distances):
    """
    Fit mixture with a reduced 4-start grid.
    Used inside bootstrap loops where we need the distribution of the
    LRT statistic, not a highly precise point estimate.  ~7x faster than
    fit_mixture for the same data.
    """
    lam0 = 1.0 / np.mean(distances)
    starts = [
        [0.1,  np.log(lam0 * 10),  np.log(lam0 * 0.2)],
        [0.1,  np.log(lam0 * 30),  np.log(lam0 * 0.2)],
        [0.05, np.log(lam0 * 10),  np.log(lam0 * 0.1)],
        [0.2,  np.log(lam0 * 20),  np.log(lam0 * 0.3)],
    ]
    best_result = None
    best_nll    = np.inf
    for x0 in starts:
        res = minimize(
            mixture_nll, x0=x0, args=(distances,), method='Nelder-Mead',
            options={'maxiter': 10000, 'xatol': 1e-6, 'fatol': 1e-6}
        )
        if res.fun < best_nll:
            best_nll    = res.fun
            best_result = res
    return best_result


def mixture_cdf(x, lam1, lam2, pi1):
    pi2 = 1.0 - pi1
    return pi1 * (1.0 - np.exp(-lam1 * x)) + pi2 * (1.0 - np.exp(-lam2 * x))


def sample_mixture(lam1, lam2, pi1, n, rng):
    """Draw n samples from the mixture via component indicator."""
    indicators = rng.random(n) < pi1
    out = np.where(
        indicators,
        rng.exponential(scale=1.0 / lam1, size=n),
        rng.exponential(scale=1.0 / lam2, size=n),
    )
    return out


def _enforce_ordering(pi1, lam1, lam2):
    """Return (pi1, lam1, lam2) with lam1 >= lam2 (high-rate first)."""
    if lam1 < lam2:
        return 1.0 - pi1, lam2, lam1
    return pi1, lam1, lam2


# ---------------------------------------------------------------------------
# Bootstrap worker functions (must be top-level for multiprocessing pickling)
# ---------------------------------------------------------------------------

def _ks_single_worker(args):
    seed, n, lam, obs_ks_stat = args
    rng = np.random.default_rng(seed)
    boot = rng.exponential(scale=1.0 / lam, size=n)
    stat, _ = kstest(boot, lambda x: single_cdf(x, lam))
    return int(stat >= obs_ks_stat)


def _ks_mixture_worker(args):
    seed, n, lam1, lam2, pi1, obs_ks_stat = args
    rng = np.random.default_rng(seed)
    boot = sample_mixture(lam1, lam2, pi1, n, rng)
    stat, _ = kstest(boot, lambda x: mixture_cdf(x, lam1, lam2, pi1))
    return int(stat >= obs_ks_stat)


def _lrt_worker(args):
    seed, n, lam_null, obs_lrt_stat = args
    rng = np.random.default_rng(seed)
    boot = rng.exponential(scale=1.0 / lam_null, size=n)
    res_s = fit_single(boot)
    res_m = fit_mixture_fast(boot)
    if res_m is None or res_m.fun >= 1e11:
        return None   # signal a skipped replicate
    lrt_boot = max(2.0 * ((-res_m.fun) - (-res_s.fun)), 0.0)
    return int(lrt_boot >= obs_lrt_stat)


# ---------------------------------------------------------------------------
# Parametric bootstrap KS p-value
# ---------------------------------------------------------------------------

def bootstrap_ks_pval_parallel(n, obs_ks_stat, worker_fn, worker_args_fn,
                                n_boot, seeds, n_workers):
    """Run a KS bootstrap in parallel. Returns p-value."""
    task_args = [worker_args_fn(s) for s in seeds]
    with multiprocessing.Pool(n_workers) as pool:
        results = pool.map(worker_fn, task_args)
    count = sum(results)
    return (count + 1) / (n_boot + 1)


# ---------------------------------------------------------------------------
# Parametric bootstrap likelihood ratio test
# ---------------------------------------------------------------------------

def bootstrap_lrt_pval_parallel(distances, obs_lrt_stat, n_boot, seeds, n_workers):
    """Run the LRT bootstrap in parallel. Returns p-value."""
    n = len(distances)
    res_null = fit_single(distances)
    lam_null = np.exp(res_null.x[0])

    task_args = [(s, n, lam_null, obs_lrt_stat) for s in seeds]
    with multiprocessing.Pool(n_workers) as pool:
        results = pool.map(_lrt_worker, task_args)

    skipped = results.count(None)
    valid   = [r for r in results if r is not None]

    if skipped > 0:
        print(f'  Warning: {skipped}/{n_boot} bootstrap LRT replicates skipped '
              f'(mixture fit failed or hit boundary).', file=sys.stderr)
    if not valid:
        print('  Warning: all bootstrap LRT replicates failed; p-value set to NaN.',
              file=sys.stderr)
        return float('nan')

    count = sum(valid)
    effective_n = len(valid)
    return (count + 1) / (effective_n + 1)


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def make_plot(distances, bins, plot_percentile,
              lam_single, lam1, lam2, pi1, output_path):
    N   = len(distances)
    pi2 = 1.0 - pi1

    xmax = np.percentile(distances, plot_percentile)
    n_truncated = int(np.sum(distances > xmax))
    if n_truncated > 0:
        print(f'  Note: plot truncated at {plot_percentile}th percentile; '
              f'{n_truncated:,} point(s) not shown.')

    counts, edges = np.histogram(distances, bins=bins, range=(0, xmax))
    bin_width   = edges[1] - edges[0]
    bin_centers = 0.5 * (edges[:-1] + edges[1:])

    x_smooth = np.linspace(0, xmax, 2000)

    single_smooth  = N * bin_width * lam_single * np.exp(-lam_single * x_smooth)
    mixture_smooth = N * bin_width * (
        pi1 * lam1 * np.exp(-lam1 * x_smooth) +
        pi2 * lam2 * np.exp(-lam2 * x_smooth)
    )
    comp1_smooth = N * bin_width * pi1 * lam1 * np.exp(-lam1 * x_smooth)
    comp2_smooth = N * bin_width * pi2 * lam2 * np.exp(-lam2 * x_smooth)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle('Exponential fits to nearest-neighbor distances', fontsize=13, y=1.01)

    for ax, yscale in zip(axes, ['linear', 'log']):
        ax.bar(bin_centers, counts, width=bin_width * 0.9,
               color='steelblue', alpha=0.45, label='data')
        ax.plot(x_smooth, single_smooth, color='orange', lw=2,
                label=f'single exp  λ={lam_single:.2e}')
        ax.plot(x_smooth, comp1_smooth, color='tomato', lw=1.5, linestyle='--',
                label=f'high-rate component  λ₁={lam1:.2e},  π₁={pi1:.3f}')
        ax.plot(x_smooth, comp2_smooth, color='seagreen', lw=1.5, linestyle='--',
                label=f'low-rate component   λ₂={lam2:.2e},  π₂={pi2:.3f}')
        ax.plot(x_smooth, mixture_smooth, color='black', lw=2,
                label='mixture fit')

        ax.set_yscale(yscale)
        ax.set_xlabel('nearest-neighbor distance (bp)')
        ax.set_ylabel('count' if yscale == 'linear' else 'count (log scale)')
        ax.set_title('linear scale' if yscale == 'linear' else 'log scale')
        ax.legend(fontsize=7.5)
        if yscale == 'log':
            ymin = max(0.5, counts[counts > 0].min() * 0.5) if np.any(counts > 0) else 0.5
            ax.set_ylim(bottom=ymin)

    info = (
        f'Fitted parameters\n'
        f'{"─"*25}\n'
        f'single:  λ = {lam_single:.4e} /bp\n'
        f'         mean gap = {1/lam_single:,.0f} bp\n'
        f'{"─"*25}\n'
        f'mixture: λ_high = {lam1:.4e} /bp\n'
        f'         λ_low  = {lam2:.4e} /bp\n'
        f'         π_high = {pi1:.4f}\n'
        f'         π_low  = {pi2:.4f}\n'
        f'         mean gap (hot)  = {1/lam1:,.0f} bp\n'
        f'         mean gap (cold) = {1/lam2:,.0f} bp'
    )
    axes[1].text(
        0.97, 0.97, info,
        transform=axes[1].transAxes,
        fontsize=7.5, verticalalignment='top', horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.85),
        fontfamily='monospace'
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f'\nPlot saved to {output_path}')
    plt.show()


# ---------------------------------------------------------------------------
# JSON results output
# ---------------------------------------------------------------------------

def write_results_json(path, distances, n_zeros_dropped, lam_single, lam1, lam2, pi1,
                       ks_stat_1, ks_pval_1,
                       ks_stat_2, ks_pval_2,
                       aic_single, bic_single,
                       aic_mixture, bic_mixture,
                       lrt_stat, lrt_pval,
                       args):
    pi2 = 1.0 - pi1
    results = {
        'script_version': __version__,
        'numpy_version':  np.__version__,
        'scipy_version':  scipy.__version__,
        'matplotlib_version': matplotlib.__version__,
        'input_file':     args.input,
        'n_distances':    len(distances),
        'n_zeros_dropped': n_zeros_dropped,
        'n_boot':         args.n_boot,
        'seed':           args.seed,
        'single_exponential': {
            'lambda':       lam_single,
            'mean_gap_bp':  1.0 / lam_single,
            'ks_statistic': ks_stat_1,
            'ks_pval_bootstrap': ks_pval_1,
            'AIC':          aic_single,
            'BIC':          bic_single,
        },
        'mixture_two_exponentials': {
            'lambda_high':      lam1,
            'lambda_low':       lam2,
            'pi_high':          pi1,
            'pi_low':           pi2,
            'mean_gap_hot_bp':  1.0 / lam1,
            'mean_gap_cold_bp': 1.0 / lam2,
            'ks_statistic':     ks_stat_2,
            'ks_pval_bootstrap': ks_pval_2,
            'AIC':              aic_mixture,
            'BIC':              bic_mixture,
        },
        'likelihood_ratio_test': {
            'LRT_statistic':    lrt_stat,
            'pval_bootstrap':   lrt_pval,
            'note': ('Single exponential is the null model. '
                     'p-value from parametric bootstrap because '
                     'null is on boundary of parameter space.')
        }
    }
    with open(path, 'w') as fh:
        json.dump(results, fh, indent=2)
    print(f'Results saved to {path}')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Fit single and mixture exponentials to nearest-neighbor distances.'
    )
    parser.add_argument('input',
                        help='Text file with one distance per line')
    parser.add_argument('--bins', type=int, default=100,
                        help='Number of histogram bins for plotting (default: 100)')
    parser.add_argument('--output', type=str, default='fit_distances.png',
                        help='Output plot filename (default: fit_distances.png)')
    parser.add_argument('--results-json', type=str, default='fit_distances.json',
                        help='Output JSON filename for fitted parameters and stats '
                             '(default: fit_distances.json)')
    parser.add_argument('--n-boot', type=int, default=999,
                        help='Bootstrap replicates for KS and LRT p-values (default: 999); '
                             'use 9999 for publication')
    parser.add_argument('--seed', type=int, default=None,
                        help='Random seed for reproducibility (default: None)')
    parser.add_argument('--plot-percentile', type=float, default=99.0,
                        help='Upper percentile of data to show in plot (default: 99.0)')
    parser.add_argument('--workers', type=int, default=multiprocessing.cpu_count(),
                        help='Parallel worker processes for bootstrap '
                             f'(default: all available CPUs)')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    args = parser.parse_args()

    # --- Library versions ---
    print(f'fit_distances.py  v{__version__}')
    print(f'  numpy {np.__version__}  |  scipy {scipy.__version__}  |  '
          f'matplotlib {matplotlib.__version__}')

    # --- RNG ---
    # Per-replicate seeds are derived from args.seed via SeedSequence later.
    if args.seed is not None:
        print(f'  random seed: {args.seed}')
    print(f'  workers: {args.workers}')

    # --- Load data ---
    try:
        distances = np.loadtxt(args.input)
    except Exception as e:
        sys.exit(f'Error reading {args.input}: {e}')

    n_raw = len(distances)
    distances = distances[distances > 0]
    n_zeros_dropped = n_raw - len(distances)
    if n_zeros_dropped > 0:
        print(f'  Dropped {n_zeros_dropped:,} zero-distance value(s).')
    if len(distances) < 10:
        sys.exit('Too few data points to fit.')

    N = len(distances)
    print(f'\nLoaded {N:,} distances.')
    print(f'  min={distances.min():.0f}  median={np.median(distances):.0f}'
          f'  max={distances.max():.0f}')

    # --- Single exponential ---
    print(f'\nFitting single exponential ...')
    res_single = fit_single(distances)
    lam_single = np.exp(res_single.x[0])
    loglik_single = -res_single.fun
    print(f'  λ = {lam_single:.4e} /bp   (mean gap = {1/lam_single:,.0f} bp)')

    # --- Mixture ---
    print(f'\nFitting mixture of two exponentials ...')
    res_mix = fit_mixture(distances)

    if res_mix is None:
        sys.exit('ERROR: Mixture fit failed entirely. '
                 'Check your data or expand starting-value grid.')
    if res_mix.fun >= 1e11:
        sys.exit('ERROR: Mixture fit converged to a boundary (pi1 outside (0,1)). '
                 'The model may not be identifiable for this dataset.')
    if not res_mix.success:
        print(f'  WARNING: optimiser flag: {res_mix.message}', file=sys.stderr)
        print(f'  Proceeding, but verify results manually.', file=sys.stderr)

    pi1_raw, log_lam1, log_lam2 = res_mix.x
    lam1_raw, lam2_raw = np.exp(log_lam1), np.exp(log_lam2)
    # Enforce lam1 >= lam2 (high-rate component first)
    pi1, lam1, lam2 = _enforce_ordering(pi1_raw, lam1_raw, lam2_raw)
    pi2 = 1.0 - pi1
    loglik_mixture = -res_mix.fun

    print(f'  λ_high = {lam1:.4e} /bp   (mean gap = {1/lam1:,.0f} bp)')
    print(f'  λ_low  = {lam2:.4e} /bp   (mean gap = {1/lam2:,.0f} bp)')
    print(f'  π_high = {pi1:.4f}')
    print(f'  π_low  = {pi2:.4f}')

    # --- AIC / BIC ---
    # Standard definitions: AIC = 2k - 2*loglik; BIC = k*ln(N) - 2*loglik
    aic_single  = 2 * 1 - 2 * loglik_single   # k=1
    aic_mixture = 2 * 3 - 2 * loglik_mixture   # k=3
    bic_single  = 1 * np.log(N) - 2 * loglik_single
    bic_mixture = 3 * np.log(N) - 2 * loglik_mixture
    delta_aic   = aic_single  - aic_mixture    # positive => mixture preferred
    delta_bic   = bic_single  - bic_mixture

    # --- LRT statistic ---
    lrt_stat = 2.0 * (loglik_mixture - loglik_single)
    if lrt_stat < 0:
        print('  WARNING: LRT statistic is negative (mixture loglik < single loglik); '
              'optimisation likely did not converge fully. Clamping to 0.',
              file=sys.stderr)
        lrt_stat = 0.0

    # --- Bootstrapped KS and LRT p-values (parallel) ---
    # Generate per-replicate seeds from the master RNG so results are
    # reproducible with --seed while each worker has independent randomness.
    ss = np.random.SeedSequence(args.seed)
    seeds = ss.spawn(3 * args.n_boot)   # 3 pools: KS-single, KS-mixture, LRT
    seeds_ks1 = [int(s.generate_state(1)[0]) for s in seeds[:args.n_boot]]
    seeds_ks2 = [int(s.generate_state(1)[0]) for s in seeds[args.n_boot:2*args.n_boot]]
    seeds_lrt = [int(s.generate_state(1)[0]) for s in seeds[2*args.n_boot:]]

    n = len(distances)
    print(f'\nRunning parametric bootstrap (n_boot={args.n_boot}, '
          f'workers={args.workers}) ...')

    print(f'  KS test: single exponential ...')
    ks_stat_1, _ = kstest(distances, lambda x: single_cdf(x, lam_single))
    ks_pval_1 = bootstrap_ks_pval_parallel(
        n, ks_stat_1,
        worker_fn=_ks_single_worker,
        worker_args_fn=lambda s: (s, n, lam_single, ks_stat_1),
        n_boot=args.n_boot, seeds=seeds_ks1, n_workers=args.workers,
    )

    print(f'  KS test: mixture ...')
    ks_stat_2, _ = kstest(distances, lambda x: mixture_cdf(x, lam1, lam2, pi1))
    ks_pval_2 = bootstrap_ks_pval_parallel(
        n, ks_stat_2,
        worker_fn=_ks_mixture_worker,
        worker_args_fn=lambda s: (s, n, lam1, lam2, pi1, ks_stat_2),
        n_boot=args.n_boot, seeds=seeds_ks2, n_workers=args.workers,
    )

    print(f'  LRT p-value (bootstrap) ...')
    lrt_pval = bootstrap_lrt_pval_parallel(
        distances, lrt_stat, args.n_boot, seeds_lrt, args.workers,
    )

    # --- Summary table ---
    print(f'\n{"─"*64}')
    print(f'{"Model":<30} {"KS stat":>10} {"bootstrap p":>12}')
    print(f'{"─"*64}')
    print(f'{"Single exponential":<30} {ks_stat_1:>10.4f} {ks_pval_1:>12.4g}')
    print(f'{"Mixture of two exponentials":<30} {ks_stat_2:>10.4f} {ks_pval_2:>12.4g}')
    print(f'{"─"*64}')

    print(f'\n{"─"*64}')
    print(f'{"Model":<30} {"AIC":>10} {"BIC":>12}')
    print(f'{"─"*64}')
    print(f'{"Single exponential":<30} {aic_single:>10.2f} {bic_single:>12.2f}')
    print(f'{"Mixture of two exponentials":<30} {aic_mixture:>10.2f} {bic_mixture:>12.2f}')
    print(f'{"─"*64}')
    preferred = 'mixture model' if delta_aic > 0 else 'single exponential'
    print(f'  ΔAIC = {abs(delta_aic):.2f} — {preferred} preferred')
    preferred = 'mixture model' if delta_bic > 0 else 'single exponential'
    print(f'  ΔBIC = {abs(delta_bic):.2f} — {preferred} preferred')
    print(f'{"─"*64}')

    print(f'\n{"─"*64}')
    print(f'Likelihood ratio test (single vs mixture)')
    print(f'  LRT statistic = {lrt_stat:.4f}')
    print(f'  bootstrap p   = {lrt_pval:.4g}')
    print(f'  (bootstrap used; null is on boundary of parameter space)')
    print(f'{"─"*64}')

    # --- Save JSON results ---
    write_results_json(
        args.results_json, distances, n_zeros_dropped,
        lam_single, lam1, lam2, pi1,
        ks_stat_1, ks_pval_1,
        ks_stat_2, ks_pval_2,
        aic_single, bic_single,
        aic_mixture, bic_mixture,
        lrt_stat, lrt_pval,
        args,
    )

    # --- Plot ---
    make_plot(distances, args.bins, args.plot_percentile,
              lam_single, lam1, lam2, pi1, args.output)


if __name__ == '__main__':
    main()
