#!/usr/bin/env python3
"""
Baum-Welch HMM for segmenting a genome into hot/cold insertion rate regions.

Model:
  - Every base is in state HOT (0) or COLD (1).
  - At each base, an insertion occurs with probability lam_hot or lam_cold.
  - The state evolves as a Markov chain with per-base transition probabilities
    p_stay_hot and p_stay_cold.
  - Each contig is treated independently.

Implementation:
  - Dense forward/backward: iterates over every base, O(genome_length).
  - Inner loops compiled with numba @njit for C-speed execution.
  - Estimated runtime: ~10-15 minutes for a 180Mb genome, 20 iterations.

Initialization:
  - STRONGLY RECOMMENDED: pass --lam_hot and --lam_cold from your mixture fit
    (fit_distances.py output). Poor initialization leads to slow convergence.

Requires: numpy, numba
    pip install numba

Input:
  - GFF file (tab-delimited). Columns: 1=contig, 4=start, 5=end (1-indexed).
    Insertion site = floor((start+end)/2), 0-indexed internally.
  - Contig lengths file: 'contig_name length' per line.

Output:
  - Fitted parameters to stdout.
  - Segmentation: 'CONTIG:START-END\tHOT/COLD' (1-indexed, inclusive).
  - With --save_posteriors: a file called 'posteriors' containing one line
    per base, tab-separated: contig, position (1-indexed), P(hot).

Usage:
    python baum_welch.py elements.gff contig_lengths.txt segmentation.txt \\
        --lam_hot 1.91e-3 --lam_cold 8.52e-5

    # Posterior decoding (default) with a stricter HOT threshold:
    python baum_welch.py elements.gff contig_lengths.txt segmentation.txt \\
        --lam_hot 1.91e-3 --lam_cold 8.52e-5 --threshold 0.9

    # Viterbi decoding (single best path) instead of posterior:
    python baum_welch.py elements.gff contig_lengths.txt segmentation.txt \\
        --lam_hot 1.91e-3 --lam_cold 8.52e-5 --viterbi

    # Save per-base posterior P(hot) to a file called "posteriors":
    python baum_welch.py elements.gff contig_lengths.txt segmentation.txt \\
        --lam_hot 1.91e-3 --lam_cold 8.52e-5 --save_posteriors

    # Re-segment from an existing "posteriors" file at a new threshold
    # (fast; skips fitting and decoding). GFF is still a required positional
    # argument but is not read in this mode.
    python baum_welch.py elements.gff contig_lengths.txt new_segmentation.txt \\
        --resegment --threshold 0.3

    python baum_welch.py --test
"""

import sys
import argparse
import numpy as np
from numba import njit


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def load_contig_lengths(path):
    lengths = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            lengths[parts[0]] = int(parts[1])
    return lengths


def load_insertions(gff_path, contig_lengths):
    insertions = {contig: [] for contig in contig_lengths}
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts  = line.strip().split('\t')
            contig = parts[0]
            start  = int(parts[3])
            end    = int(parts[4])
            idx    = (start + end) // 2 - 1  # midpoint, 0-indexed
            if contig not in insertions:
                print(f'Warning: {contig} not in lengths file, skipping.')
                continue
            insertions[contig].append(idx)
    for contig in insertions:
        insertions[contig] = sorted(insertions[contig])
    return insertions


def build_obs(ins_positions, contig_length):
    """Build binary observation array: 1 at each insertion, 0 elsewhere."""
    obs = np.zeros(contig_length, dtype=np.float64)
    for pos in ins_positions:
        if 0 <= pos < contig_length:
            obs[pos] = 1.0
    return obs


# ---------------------------------------------------------------------------
# Forward pass (numba-compiled)
#
# alpha[t, s] = P(obs[0..t], state[t]=s), rescaled at each step.
# log_norm    = sum of log scale factors = log P(all observations).
# ---------------------------------------------------------------------------

@njit
def forward(obs, pi, A, lam):
    n     = len(obs)
    alpha = np.zeros((n, 2))

    # Initialise at t=0
    for s in range(2):
        if obs[0] == 1.0:
            alpha[0, s] = pi[s] * lam[s]
        else:
            alpha[0, s] = pi[s] * (1.0 - lam[s])
    sc = alpha[0, 0] + alpha[0, 1]
    alpha[0, 0] /= sc
    alpha[0, 1] /= sc
    log_norm = np.log(sc)

    # Recursion
    for t in range(1, n):
        for s in range(2):
            val = 0.0
            for prev in range(2):
                val += alpha[t-1, prev] * A[prev, s]
            if obs[t] == 1.0:
                val *= lam[s]
            else:
                val *= (1.0 - lam[s])
            alpha[t, s] = val
        sc = alpha[t, 0] + alpha[t, 1]
        if sc > 0.0:
            alpha[t, 0] /= sc
            alpha[t, 1] /= sc
            log_norm += np.log(sc)

    return alpha, log_norm


# ---------------------------------------------------------------------------
# Backward pass (numba-compiled)
#
# beta[t, s] = P(obs[t+1..n-1] | state[t]=s), rescaled at each step.
# ---------------------------------------------------------------------------

@njit
def backward(obs, A, lam):
    n    = len(obs)
    beta = np.ones((n, 2))

    for t in range(n - 2, -1, -1):
        for s in range(2):
            val = 0.0
            for nxt in range(2):
                if obs[t+1] == 1.0:
                    emit = lam[nxt]
                else:
                    emit = 1.0 - lam[nxt]
                val += A[s, nxt] * emit * beta[t+1, nxt]
            beta[t, s] = val
        sc = beta[t, 0] + beta[t, 1]
        if sc > 0.0:
            beta[t, 0] /= sc
            beta[t, 1] /= sc

    return beta


# ---------------------------------------------------------------------------
# Sufficient statistics (numba-compiled)
#
# gamma_emit[s]  = expected insertions in state s
# gamma_bases[s] = expected bases in state s
# xi_sum[i,j]    = expected transitions from state i to state j
# gamma_init[s]  = posterior state probability at t=0 (for pi update)
#
# alpha and beta are independently rescaled, but their product gives
# correct RELATIVE posteriors (scale factors cancel in the ratio).
# ---------------------------------------------------------------------------

@njit
def compute_stats(obs, alpha, beta, A, lam):
    n           = len(obs)
    gamma_emit  = np.zeros(2)
    gamma_bases = np.zeros(2)
    xi_sum      = np.zeros((2, 2))

    for t in range(n):
        # Gamma at time t
        raw = alpha[t, 0] * beta[t, 0] + alpha[t, 1] * beta[t, 1]
        if raw > 0.0:
            g0 = alpha[t, 0] * beta[t, 0] / raw
            g1 = alpha[t, 1] * beta[t, 1] / raw
        else:
            g0 = 0.5
            g1 = 0.5

        gamma_bases[0] += g0
        gamma_bases[1] += g1
        if obs[t] == 1.0:
            gamma_emit[0] += g0
            gamma_emit[1] += g1

        # Xi at time t (transition to t+1)
        if t < n - 1:
            total  = 0.0
            xi_raw = np.zeros((2, 2))
            for i in range(2):
                for j in range(2):
                    emit = lam[j] if obs[t+1] == 1.0 else 1.0 - lam[j]
                    xi_raw[i, j] = alpha[t, i] * A[i, j] * emit * beta[t+1, j]
                    total += xi_raw[i, j]
            if total > 0.0:
                for i in range(2):
                    for j in range(2):
                        xi_sum[i, j] += xi_raw[i, j] / total

    # Initial state
    raw        = alpha[0, 0] * beta[0, 0] + alpha[0, 1] * beta[0, 1]
    gamma_init = np.zeros(2)
    if raw > 0.0:
        gamma_init[0] = alpha[0, 0] * beta[0, 0] / raw
        gamma_init[1] = alpha[0, 1] * beta[0, 1] / raw
    else:
        gamma_init[0] = 0.5
        gamma_init[1] = 0.5

    return gamma_emit, gamma_bases, xi_sum, gamma_init


# ---------------------------------------------------------------------------
# Viterbi (numba-compiled)
# ---------------------------------------------------------------------------

@njit
def viterbi(obs, pi, A, lam):
    n     = len(obs)
    delta = np.zeros((n, 2))
    psi   = np.zeros((n, 2), dtype=np.int64)

    for s in range(2):
        emit         = np.log(lam[s] + 1e-300) if obs[0] == 1.0 \
                       else np.log(1.0 - lam[s] + 1e-300)
        delta[0, s]  = np.log(pi[s] + 1e-300) + emit

    for t in range(1, n):
        for s in range(2):
            best_val  = -1e300
            best_prev = 0
            for prev in range(2):
                val = delta[t-1, prev] + np.log(A[prev, s] + 1e-300)
                if val > best_val:
                    best_val  = val
                    best_prev = prev
            emit         = np.log(lam[s] + 1e-300) if obs[t] == 1.0 \
                           else np.log(1.0 - lam[s] + 1e-300)
            delta[t, s]  = best_val + emit
            psi[t, s]    = best_prev

    states     = np.zeros(n, dtype=np.int64)
    states[-1] = 0 if delta[n-1, 0] >= delta[n-1, 1] else 1
    for t in range(n - 2, -1, -1):
        states[t] = psi[t+1, states[t+1]]

    return states


# ---------------------------------------------------------------------------
# Posterior decoding (numba-compiled)
#
# Computes gamma[t, s] = P(state[t]=s | all observations) at every position,
# then calls state HOT if gamma[t, 0] > threshold (default 0.5).
#
# More sensitive than Viterbi: integrates over all paths rather than
# committing to the single best path. A cluster of insertions that is not
# on the globally best Viterbi path may still have high posterior probability.
# ---------------------------------------------------------------------------

@njit
def posterior_decode(alpha, beta, threshold):
    n      = len(alpha)
    states = np.zeros(n, dtype=np.int64)
    for t in range(n):
        raw = alpha[t, 0] * beta[t, 0] + alpha[t, 1] * beta[t, 1]
        if raw > 0.0:
            p_hot = alpha[t, 0] * beta[t, 0] / raw
        else:
            p_hot = 0.5
        states[t] = 0 if p_hot > threshold else 1
    return states


@njit
def compute_p_hot(alpha, beta):
    """Return per-base P(state = hot | all observations)."""
    n     = len(alpha)
    p_hot = np.zeros(n)
    for t in range(n):
        raw = alpha[t, 0] * beta[t, 0] + alpha[t, 1] * beta[t, 1]
        if raw > 0.0:
            p_hot[t] = alpha[t, 0] * beta[t, 0] / raw
        else:
            p_hot[t] = 0.5
    return p_hot


# ---------------------------------------------------------------------------
# Warmup: compile numba functions on a tiny array
# ---------------------------------------------------------------------------

def warmup():
    print('Compiling numba functions (one-time, ~10s) ...')
    obs_tiny = np.zeros(100, dtype=np.float64)
    obs_tiny[10] = 1.0; obs_tiny[50] = 1.0
    pi  = np.array([0.5, 0.5])
    A   = np.array([[0.9999, 0.0001],[0.0001, 0.9999]])
    lam = np.array([1e-3, 1e-5])
    a, _ = forward(obs_tiny, pi, A, lam)
    b    = backward(obs_tiny, A, lam)
    compute_stats(obs_tiny, a, b, A, lam)
    viterbi(obs_tiny, pi, A, lam)
    posterior_decode(a, b, 0.5)
    compute_p_hot(a, b)
    print('Done.\n')


# ---------------------------------------------------------------------------
# Baum-Welch
# ---------------------------------------------------------------------------

def baum_welch(insertions, contig_lengths,
               lam_hot=1e-3, lam_cold=1e-5,
               p_stay_hot=0.9999, p_stay_cold=0.9999,
               max_iter=100, tol=1e-6):
    import time

    # pi is fixed to the stationary distribution implied by the current
    # transition matrix (recomputed each iteration below). This avoids the
    # feedback loop where EM-estimated pi from position 0 of each contig
    # drifts toward whatever state early iterations happened to assign to
    # contig starts. Contig breaks are arbitrary, so there is no reason
    # to learn a special initial distribution.
    pi = np.array([
        (1 - p_stay_cold) / ((1 - p_stay_hot) + (1 - p_stay_cold)),
        (1 - p_stay_hot)  / ((1 - p_stay_hot) + (1 - p_stay_cold)),
    ])
    prev_log_lik = -np.inf

    print(f'Initial parameters:')
    print(f'  lam_hot     = {lam_hot:.4e}')
    print(f'  lam_cold    = {lam_cold:.4e}')
    print(f'  p_stay_hot  = {p_stay_hot:.6f}')
    print(f'  p_stay_cold = {p_stay_cold:.6f}')
    print()

    for iteration in range(max_iter):
        A   = np.array([[p_stay_hot,       1 - p_stay_hot ],
                        [1 - p_stay_cold,  p_stay_cold    ]])
        lam = np.array([lam_hot, lam_cold])

        total_gamma_emit  = np.zeros(2)
        total_gamma_bases = np.zeros(2)
        total_xi_sum      = np.zeros((2, 2))
        total_gamma_init  = np.zeros(2)
        total_log_lik     = 0.0
        t0 = time.time()

        for contig, ins_positions in insertions.items():
            L   = contig_lengths[contig]
            obs = build_obs(ins_positions, L)

            alpha, log_norm = forward(obs, pi, A, lam)
            total_log_lik  += log_norm

            beta = backward(obs, A, lam)

            ge, gb, xi, gi = compute_stats(obs, alpha, beta, A, lam)

            total_gamma_emit  += ge
            total_gamma_bases += gb
            total_xi_sum      += xi
            total_gamma_init  += gi

        # M step
        for i in range(2):
            row = np.sum(total_xi_sum[i, :]) + 1e-300
            A[i, :] = total_xi_sum[i, :] / row

        p_stay_hot  = float(np.clip(A[0, 0], 1e-8, 1 - 1e-8))
        p_stay_cold = float(np.clip(A[1, 1], 1e-8, 1 - 1e-8))
        lam_hot     = float(np.clip(
            total_gamma_emit[0] / (total_gamma_bases[0] + 1e-300), 1e-8, 0.5))
        lam_cold    = float(np.clip(
            total_gamma_emit[1] / (total_gamma_bases[1] + 1e-300), 1e-8, 0.5))

        # pi recomputed from the updated transition matrix (stationary
        # distribution), not estimated from posterior at position 0.
        pi = np.array([
            (1 - p_stay_cold) / ((1 - p_stay_hot) + (1 - p_stay_cold)),
            (1 - p_stay_hot)  / ((1 - p_stay_hot) + (1 - p_stay_cold)),
        ])

        elapsed = time.time() - t0
        print(f'Iteration {iteration + 1:3d} ({elapsed:.1f}s): '
              f'log-lik={total_log_lik:.2f}  '
              f'lam_hot={lam_hot:.4e}  lam_cold={lam_cold:.4e}  '
              f'p_stay_hot={p_stay_hot:.8f}  p_stay_cold={p_stay_cold:.8f}')

        if abs(total_log_lik - prev_log_lik) < tol:
            print(f'Converged after {iteration + 1} iterations.')
            break

        prev_log_lik = total_log_lik

    return lam_hot, lam_cold, p_stay_hot, p_stay_cold, pi, A


# ---------------------------------------------------------------------------
# Convert Viterbi state sequence to genomic regions
# ---------------------------------------------------------------------------

def states_to_regions(states, contig, contig_length):
    regions = []
    current = states[0]
    start   = 0  # 0-indexed

    for t in range(1, len(states)):
        if states[t] != current:
            label = 'HOT' if current == 0 else 'COLD'
            regions.append(f'{contig}:{start + 1}-{t}\t{label}')
            current = states[t]
            start   = t

    label = 'HOT' if current == 0 else 'COLD'
    regions.append(f'{contig}:{start + 1}-{contig_length}\t{label}')
    return regions


# ---------------------------------------------------------------------------
# Simulation for testing
# ---------------------------------------------------------------------------

def simulate(L, lam_hot, lam_cold, p_stay_hot, p_stay_cold, seed=42):
    rng   = np.random.default_rng(seed)
    state = 0
    ins   = []
    for t in range(L):
        lam = lam_hot if state == 0 else lam_cold
        if rng.random() < lam:
            ins.append(t)
        state = (0 if rng.random() < p_stay_hot  else 1) if state == 0 \
                else (1 if rng.random() < p_stay_cold else 0)
    return ins


# ---------------------------------------------------------------------------
# Re-segment from a saved posteriors file
#
# Reads the 'posteriors' file (contig, position, p_hot per line), bins
# posteriors by contig, and writes a new segmentation file at a chosen
# threshold. Avoids rerunning forward-backward.
# ---------------------------------------------------------------------------

def resegment_from_posteriors(posteriors_path, contig_lengths, threshold, output_path):
    import time
    print(f'Reading posteriors from {posteriors_path} ...')
    t0 = time.time()

    # Accumulate p_hot arrays per contig. We allocate the full array up front
    # and fill positions as we read them, so we don't care about line order.
    arrays = {c: np.full(L, np.nan) for c, L in contig_lengths.items()}

    n_lines = 0
    with open(posteriors_path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 3:
                continue
            contig = parts[0]
            pos    = int(parts[1]) - 1  # 1-indexed in file -> 0-indexed
            p_hot  = float(parts[2])
            if contig in arrays and 0 <= pos < len(arrays[contig]):
                arrays[contig][pos] = p_hot
                n_lines += 1

    print(f'  Read {n_lines:,} posterior values in {time.time()-t0:.1f}s')

    # Check for any unfilled positions (NaN). If so, warn.
    missing = {c: int(np.isnan(arr).sum()) for c, arr in arrays.items()}
    total_missing = sum(missing.values())
    if total_missing > 0:
        print(f'  Warning: {total_missing:,} positions across {sum(1 for v in missing.values() if v > 0)} '
              f'contigs had no posterior value (treated as cold).')

    print(f'Thresholding at {threshold} and writing segmentation ...')
    with open(output_path, 'w') as out:
        for contig, L in contig_lengths.items():
            p_hot_arr = arrays[contig]
            # Threshold: hot (state 0) if p_hot > threshold, else cold (state 1).
            # NaN positions are treated as cold (p_hot > threshold is False for NaN).
            states = np.where(p_hot_arr > threshold, 0, 1).astype(np.int64)
            for line in states_to_regions(states, contig, L):
                out.write(line + '\n')

    print(f'Segmentation written to {output_path}')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) == 2 and sys.argv[1] == '--test':
        print('Running simulation test ...\n')
        L             = 5000000
        true_lam_hot  = 1e-3
        true_lam_cold = 1e-5
        true_psh      = 0.9990
        true_psc      = 0.9995
        print(f'True parameters:')
        print(f'  lam_hot     = {true_lam_hot:.4e}')
        print(f'  lam_cold    = {true_lam_cold:.4e}')
        print(f'  p_stay_hot  = {true_psh:.6f}')
        print(f'  p_stay_cold = {true_psc:.6f}')
        print()
        ins = simulate(L, true_lam_hot, true_lam_cold, true_psh, true_psc)
        print(f'Simulated {len(ins)} insertions on {L:,} bp\n')
        warmup()
        insertions     = {'test': ins}
        contig_lengths = {'test': L}
        lh, lc, psh, psc, pi, A = baum_welch(
            insertions, contig_lengths,
            lam_hot=true_lam_hot, lam_cold=true_lam_cold,
            p_stay_hot=true_psh, p_stay_cold=true_psc)
        print(f'\nRecovered parameters:')
        print(f'  lam_hot     = {lh:.4e}  (true: {true_lam_hot:.4e})')
        print(f'  lam_cold    = {lc:.4e}  (true: {true_lam_cold:.4e})')
        print(f'  p_stay_hot  = {psh:.6f}  (true: {true_psh:.6f})')
        print(f'  p_stay_cold = {psc:.6f}  (true: {true_psc:.6f})')
        return

    parser = argparse.ArgumentParser(
        description='Baum-Welch HMM segmentation of insertion rate regions.')
    parser.add_argument('gff',     help='GFF file of elements')
    parser.add_argument('lengths', help='Contig lengths file')
    parser.add_argument('output',  help='Output segmentation file')
    parser.add_argument('--lam_hot',     type=float, default=1e-3,
        help='Initial lam_hot (default: 1e-3). Use value from fit_distances.py.')
    parser.add_argument('--lam_cold',    type=float, default=1e-5,
        help='Initial lam_cold (default: 1e-5). Use value from fit_distances.py.')
    parser.add_argument('--p_stay_hot',  type=float, default=0.9999,
        help='Initial p_stay_hot (default: 0.9999)')
    parser.add_argument('--p_stay_cold', type=float, default=0.9999,
        help='Initial p_stay_cold (default: 0.9999)')
    parser.add_argument('--max_iter',    type=int,   default=100,
        help='Maximum Baum-Welch iterations (default: 100)')
    parser.add_argument('--tol',         type=float, default=1e-6,
        help='Convergence tolerance on log-likelihood (default: 1e-6)')
    parser.add_argument('--threshold',   type=float, default=0.5,
        help='Posterior probability threshold for HOT call (default: 0.5)')
    parser.add_argument('--viterbi',     action='store_true',
        help='Use Viterbi decoding instead of posterior decoding')
    parser.add_argument('--save_posteriors', action='store_true',
        help='Write per-base posterior P(hot) to a file called "posteriors" '
             '(tab-delimited: contig, position (1-indexed), p_hot). Incompatible '
             'with --viterbi. File is large (~5-8 GB for a 180Mb genome).')
    parser.add_argument('--resegment', action='store_true',
        help='Skip fitting and decoding; read an existing "posteriors" file and '
             'produce a new segmentation using --threshold. Fast (seconds). '
             'The GFF and initial parameters are ignored in this mode, but the '
             'contig lengths file is still required and must match the one used '
             'to produce the posteriors.')
    args = parser.parse_args()

    # Argument validation
    if args.save_posteriors and args.viterbi:
        parser.error('--save_posteriors is incompatible with --viterbi '
                     '(Viterbi does not compute posteriors).')

    # Resegment-from-file mode: skip fitting and decoding, just threshold
    # an existing posteriors file.
    if args.resegment:
        print('Loading contig lengths ...')
        contig_lengths = load_contig_lengths(args.lengths)
        print(f'  {len(contig_lengths)} contigs, total {sum(contig_lengths.values()):,} bp')
        resegment_from_posteriors('posteriors', contig_lengths,
                                  args.threshold, args.output)
        return

    print('Loading contig lengths ...')
    contig_lengths = load_contig_lengths(args.lengths)
    print(f'  {len(contig_lengths)} contigs, total {sum(contig_lengths.values()):,} bp')

    print('Loading GFF ...')
    insertions = load_insertions(args.gff, contig_lengths)
    print(f'  {sum(len(v) for v in insertions.values()):,} insertions\n')

    warmup()

    print('Running Baum-Welch ...')
    lh, lc, psh, psc, pi, A = baum_welch(
        insertions, contig_lengths,
        lam_hot=args.lam_hot, lam_cold=args.lam_cold,
        p_stay_hot=args.p_stay_hot, p_stay_cold=args.p_stay_cold,
        max_iter=args.max_iter, tol=args.tol)

    # Stationary distribution from transition matrix (more reliable than EM pi)
    switch_to_cold    = 1 - psh
    switch_to_hot     = 1 - psc
    pi_hot_stationary = switch_to_hot  / (switch_to_hot + switch_to_cold)
    pi_cold_stationary= switch_to_cold / (switch_to_hot + switch_to_cold)

    print(f'\nFitted parameters:')
    print(f'  lam_hot     = {lh:.4e}  (mean gap = {1/lh:,.0f} bp)')
    print(f'  lam_cold    = {lc:.4e}  (mean gap = {1/lc:,.0f} bp)')
    print(f'  p_stay_hot  = {psh:.8f}  (mean hot  region = {1/(1-psh):,.0f} bp)')
    print(f'  p_stay_cold = {psc:.8f}  (mean cold region = {1/(1-psc):,.0f} bp)')
    print(f'  pi_hot      = {pi_hot_stationary:.4f}  (from transition matrix)')
    print(f'  pi_cold     = {pi_cold_stationary:.4f}  (from transition matrix)')

    decoding = 'Viterbi' if args.viterbi else f'posterior (threshold={args.threshold})'
    print(f'\nDecoding with {decoding} and writing segmentation ...')
    if args.save_posteriors:
        print('Also writing per-base posteriors to "posteriors" '
              '(this will be large).')
    lam = np.array([lh, lc])
    # Use stationary pi for decoding (treat each contig start as a random
    # draw from the long-run state distribution, not as a distinguished
    # position with its own learned prior).
    pi = np.array([pi_hot_stationary, pi_cold_stationary])
    post_fh = open('posteriors', 'w') if args.save_posteriors else None
    try:
        with open(args.output, 'w') as out:
            for contig, ins_positions in insertions.items():
                L   = contig_lengths[contig]
                obs = build_obs(ins_positions, L)
                if args.viterbi:
                    states = viterbi(obs, pi, A, lam)
                else:
                    alpha, _ = forward(obs, pi, A, lam)
                    beta     = backward(obs, A, lam)
                    states   = posterior_decode(alpha, beta, args.threshold)
                    if post_fh is not None:
                        p_hot = compute_p_hot(alpha, beta)
                        # Write every base: contig<TAB>pos(1-indexed)<TAB>p_hot
                        for t in range(L):
                            post_fh.write(f'{contig}\t{t + 1}\t{p_hot[t]:.6f}\n')
                for line in states_to_regions(states, contig, L):
                    out.write(line + '\n')
    finally:
        if post_fh is not None:
            post_fh.close()

    print(f'Segmentation written to {args.output}')
    if args.save_posteriors:
        print('Per-base posteriors written to "posteriors".')


if __name__ == '__main__':
    main()
