#!/usr/bin/env python3.10
"""Monte Carlo of the reduced-scale SNOW-Axn analogue.

Draws N i.i.d. initial states (uniform 64-bit LFSR + 24-bit FSM), clocks the
cipher T = 8 times and measures, for each of the 7 nonzero masks of the
3-dimensional span V_3 and the clocks t = 1..6,

  (i)   single-clock correlation   tau_emp = E[(-1)^{e(t)}]      (t = 0 too),
  (ii)  folded correlation         C2(d)    = E[(-1)^{e(t) ^ e(t+d)}],
  (iii) the different-clock independence test C2(d) vs tau_emp(t)tau_emp(t+d),
  (iv)  the empirical capacity     SEI_2 = sum_w C2(w)^2  (d = 1),
  (v)   the t = 0 sanity check (independent uniform FSM registers).

Standard errors are batch means: the N samples are split into B independent
batches, the quantity is recomputed on each batch and SE = std(batch
estimates)/sqrt(B).  This is robust against the correlation between clocks
inside one sample.

The masks (3 linearly independent triples of |tau| = 2^-9, and their 7 nonzero
span combinations) are loaded from masks.json when present, otherwise the
search is recomputed.

Run (from any directory):
    python3.10 monte_carlo.py --n 24            # N = 2^24
    python3.10 monte_carlo.py --n 22 --json out.json
"""

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_tables as bt
import search_masks as sm

HERE = os.path.dirname(os.path.abspath(__file__))
T_CLOCKS = 8          # clock 0..7, e(t) measured for t = 1..6
E_MIN, E_MAX = 1, 6   # clocks on which e(t) is measured


def load_basis(path=os.path.join(HERE, "masks.json"), v=3):
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            payload = json.load(f)
        basis = [(r["a"], r["b"], r["c"]) for r in payload["basis"]]
        if len(basis) >= v:
            span = [(r["w"], r["a"], r["b"], r["c"], r["members"])
                    for r in payload["span"]]
            return basis[:v], span
    T = bt.get_tables()
    top = sm.find_top(T, 40)
    chosen = sm.independent_basis(top, v)
    basis = [(a, b, c) for (_, a, b, c) in chosen]
    span = sm.span_masks(chosen)
    return basis, span


def combo_masks(basis):
    """Return list (w, a, b, c) of the 7 nonzero span combinations."""
    out = []
    for w in range(1, 8):
        a = b = c = 0
        for i, (ai, bi, ci) in enumerate(basis):
            if (w >> i) & 1:
                a ^= ai
                b ^= bi
                c ^= ci
        out.append((w, a, b, c))
    return out


def run_mc(n, seed, batches, chunk, basis, span):
    T = bt.get_tables()
    rng = np.random.default_rng(seed)
    per_batch = n // batches
    assert per_batch * batches == n
    ncombos = 7
    nclk = E_MAX + 1                       # e(t) for t = 0..6
    combos = combo_masks(basis)            # [(w,a,b,c)] in w = 1..7 order
    base_masks = [(a, b, c) for (a, b, c) in basis]       # the 3 basis triples
    singles = np.zeros((batches, ncombos, nclk), dtype=np.int64)
    paircnt = np.zeros((batches, ncombos, 3, E_MAX + 1), dtype=np.int64)
    for b in range(batches):
        left = per_batch
        while left:
            m = min(chunk, left)
            left -= m
            u = [rng.integers(0, 256, size=m, dtype=np.uint8) for _ in range(8)]
            r1 = rng.integers(0, 256, size=m, dtype=np.uint8)
            r2 = rng.integers(0, 256, size=m, dtype=np.uint8)
            r3 = rng.integers(0, 256, size=m, dtype=np.uint8)
            # H(t), t = 0..7
            H = [u[0]]
            for step in range(1, T_CLOCKS):
                u = T.lfsr_advance(u)
                H.append(u[0])
            # register states at t = 0..7
            states = [(r1, r2, r3)]
            for step in range(1, T_CLOCKS):
                r1, r2, r3 = T.fsm_advance(r1, r2, r3, H[step])
                states.append((r1, r2, r3))
            # e_w(t) for t = 0..6, w = 1..7 (combos = XOR of the base triples)
            E = np.zeros((nclk, ncombos, m), dtype=np.uint8)
            for t in range(nclk):
                A, B, C = T.block_vectors(*states[t])
                base = [T.parity[(a & A) ^ (b & B) ^ (c & C)]
                        for (a, b, c) in base_masks]
                for wi in range(ncombos):
                    w = wi + 1
                    arr = None
                    for i in range(3):
                        if (w >> i) & 1:
                            arr = base[i] if arr is None else (arr ^ base[i])
                    E[t, wi] = arr
            for t in range(nclk):
                singles[b, :, t] += E[t].sum(axis=1, dtype=np.int64)
            for d in (1, 2, 3):
                for t in range(E_MIN, E_MAX - d + 1):
                    paircnt[b, :, d - 1, t] += (E[t] ^ E[t + d]).sum(
                        axis=1, dtype=np.int64)
    # -- batch estimates ------------------------------------------------------
    single_est = 1.0 - 2.0 * singles / float(per_batch)        # [B,7,7]
    folded_est = 1.0 - 2.0 * paircnt / float(per_batch)        # [B,7,3,7]
    return single_est, folded_est, per_batch


def se(x, axis=0):
    """SE of the mean over independent batches."""
    return x.std(axis=axis, ddof=1) / np.sqrt(x.shape[axis])


def summarise(single_est, folded_est, basis, span, n, batches, runtime):
    T = bt.get_tables()
    combos = combo_masks(basis)
    ncombos = 7
    out = {"n": n, "batches": batches, "runtime_s": runtime,
           "expected_se": 1.0 / np.sqrt(n), "combos": [], "span": [],
           "basis": []}

    # predicted tau for every combo
    tau_pred = [bt.tau_unit(a, b, c, T) / 16777216.0 for (_, a, b, c) in combos]
    for i, (w, a, b, c) in enumerate(combos):
        out["span"].append({"w": w, "a": a, "b": b, "c": c,
                            "cprime": bt.sigma_inv(c),
                            "tau_pred": tau_pred[i]})

    print("#" * 78)
    print("# Reduced-scale SNOW-Axn analogue -- Monte Carlo")
    print("# N = 2^%.3f (%d samples), batches = %d, runtime = %.1f s"
          % (np.log2(n), n, batches, runtime))
    print("# sigma = rotl3, LFSR 64 bit + FSM 24 bit uniform initial state,")
    print("# T = %d clocks, e(t) measured at t = %d..%d" % (T_CLOCKS, E_MIN, E_MAX))
    print("# SE of a single mean from N samples: 1/sqrt(N) = %.3e"
          % (1.0 / np.sqrt(n)))
    print("#" * 78)

    # -- (i) single-clock correlations --------------------------------------
    print()
    print("# (i) single-clock correlation  tau_emp = E[(-1)^{e(t)}]")
    print("#     a     b     c   | t=0 (indep.)      t=1..6 mean         SE        | tau_pred")
    for i, (w, a, b, c) in enumerate(combos):
        s0 = single_est[:, i, 0]
        s0m, s0e = s0.mean(), se(s0)
        st = single_est[:, i, E_MIN:E_MAX + 1].mean(axis=1)
        stm, ste = st.mean(), se(st)
        print("  0x%02x  0x%02x  0x%02x | %+0.9f  %+0.9f  %+0.9f  | %+0.9f"
              % (a, b, c, s0m, stm, ste, tau_pred[i]))
        out["combos"].append({
            "w": w, "a": a, "b": b, "c": c, "tau_pred": tau_pred[i],
            "single_t0": s0m, "single_t0_se": s0e,
            "single_t1_6": stm, "single_t1_6_se": ste,
            "single_per_t": [float(v) for v in
                             single_est[:, i, E_MIN:E_MAX + 1].mean(axis=0)],
        })

    # -- (ii) folded correlations -------------------------------------------
    print()
    print("# (ii) folded check correlation  C2(d) = E[(-1)^{e(t) ^ e(t+d)}]")
    print("#      prediction under independent two-clock noise: tau_pred^2")
    print("#   combo        d | empirical C2        SE          tau_pred^2    diff/SE")
    folded = {}
    for i, (w, a, b, c) in enumerate(combos):
        row = []
        for d in (1, 2, 3):
            est = folded_est[:, i, d - 1, E_MIN:E_MAX - d + 1].mean(axis=1)
            m, e = est.mean(), se(est)
            pred = tau_pred[i] ** 2
            diff = m - pred
            row.append({"d": d, "C2": m, "se": e, "tau2_pred": pred,
                        "diff_over_se": (diff / e if e else 0.0),
                        "per_t": [float(v) for v in
                                  folded_est[:, i, d - 1, E_MIN:E_MAX - d + 1].mean(axis=0)]})
            print("  0x%02x 0x%02x 0x%02x  %d | %+0.9f  %+0.9f  %+0.9e  %+10.2f"
                  % (a, b, c, d, m, e, pred, (diff / e if e else 0.0)))
        folded[w] = row
        out["combos"][i]["folded"] = row

    # -- (iii) different-clock independence ---------------------------------
    print()
    print("# (iii) different-clock independence: C2(d) vs tau_emp(t)*tau_emp(t+d)")
    print("#       (difference and significance; df/SE > ~3 means a real")
    print("#        correlation between the two clocks' approximation errors)")
    print("#   combo        d | C2(d)             tau(t)*tau(t+d)     difference      SE      diff/SE")
    for i, (w, a, b, c) in enumerate(combos):
        for d in (1, 2, 3):
            ts = list(range(E_MIN, E_MAX - d + 1))
            est = folded_est[:, i, d - 1, E_MIN:E_MAX - d + 1].mean(axis=1)
            prod = np.ones(folded_est.shape[0])
            for t in ts:
                prod *= single_est[:, i, t]
            diff = est - prod
            m, e, dm, de = est.mean(), se(est), diff.mean(), se(diff)
            print("  0x%02x 0x%02x 0x%02x  %d | %+0.9f  %+0.9f  %+0.9e  %+0.9e  %+8.2f"
                  % (a, b, c, d, m, prod.mean(), dm, de, (dm / de if de else 0.0)))
            out["combos"][i].setdefault("independence", []).append(
                {"d": d, "C2": m, "C2_se": e, "single_product": float(prod.mean()),
                 "difference": dm, "difference_se": de,
                 "diff_over_se": (dm / de if de else 0.0)})

    # -- (iv) capacity ------------------------------------------------------
    print()
    print("# (iv) empirical capacity SEI_2 = sum_w C2(w)^2  at d = 1")
    c2_batch = np.stack([folded_est[:, i, 0, E_MIN:E_MAX].mean(axis=1)
                         for i in range(ncombos)], axis=1)      # [B,7]
    sei_batch = (c2_batch ** 2).sum(axis=1)
    sei_mean, sei_se = sei_batch.mean(), se(sei_batch)
    pred = sum(t ** 4 for t in tau_pred)
    print("#  empirical  SEI_2 = %.9e +/- %.3e   (log2 = %.6f +/- %.4f)"
          % (sei_mean, sei_se, np.log2(sei_mean), sei_se / sei_mean / np.log(2)))
    print("#  predicted  SEI_2 = %.9e            (log2 = %.6f)"
          % (pred, np.log2(pred)))
    print("#  empirical / predicted = %.3f  (log2 difference = %+.3f bits)"
          % (sei_mean / pred, np.log2(sei_mean / pred)))
    print("#  per-combo C2(w)^2:")
    for i, (w, a, b, c) in enumerate(combos):
        print("#   w=%d (0x%02x,0x%02x,0x%02x): C2=%.6e  C2^2=%.6e  tau_pred^2=%.6e"
              % (w, a, b, c, c2_batch[:, i].mean(), (c2_batch[:, i].mean()) ** 2,
                 tau_pred[i] ** 2))
    out["capacity"] = {"empirical_SEI2": sei_mean, "empirical_SEI2_se": sei_se,
                       "predicted_SEI2": pred,
                       "log2_empirical": float(np.log2(sei_mean)),
                       "log2_predicted": float(np.log2(pred)),
                       "log2_ratio": float(np.log2(sei_mean / pred))}

    # -- summary of t = 0 vs t > 0 -----------------------------------------
    print()
    print("# summary")
    print("#  t=0 (independent uniform FSM): mean |tau_emp - tau_pred| over combos =")
    dev0 = np.mean([abs(out["combos"][i]["single_t0"] - tau_pred[i])
                    for i in range(ncombos)])
    se0 = 1.0 / np.sqrt(n)
    print("#     %.3e  (expected MC deviation ~ %.3e = 1/sqrt(N))" % (dev0, se0))
    print("#  t=1..6: mean |tau_emp - tau_pred| over combos =")
    devt = np.mean([abs(out["combos"][i]["single_t1_6"] - tau_pred[i])
                    for i in range(ncombos)])
    print("#     %.3e" % devt)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=22,
                    help="log2 of the sample count (22 -> 2^22, 24 -> 2^24)")
    ap.add_argument("--seed", type=int, default=20261121)
    ap.add_argument("--batches", type=int, default=32)
    ap.add_argument("--chunk", type=int, default=1 << 20)
    ap.add_argument("--json", default=None)
    args = ap.parse_args()
    n = 1 << args.n
    basis, span = load_basis()
    t0 = time.time()
    single_est, folded_est, per_batch = run_mc(n, args.seed, args.batches,
                                              args.chunk, basis, span)
    runtime = time.time() - t0
    out = summarise(single_est, folded_est, basis, span, n, args.batches, runtime)
    if args.json:
        with open(args.json, "w", encoding="utf-8") as f:
            json.dump(out, f, indent=1)
        print("# wrote", args.json)


if __name__ == "__main__":
    main()
