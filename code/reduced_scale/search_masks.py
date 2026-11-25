#!/usr/bin/env python3.10
"""Search the linear masks of the reduced-scale SNOW-Axn analogue.

We maximise |tau(a,b,c)| over all triples of nonzero byte masks.  With
c' = sigma^-1(c) and LAT[a][b] = sum_x (-1)^{a.x ^ b.S(x)},

    tau(a,b,c) = LAT[b][c] * LAT[a][b^c'] * LAT[a][c'] / 2^24 .

(The three factors are the correlations of the three one-register blocks:
e_R1(x) = b.x ^ c.S(x), e_R2(x) = a.S^-1(x) ^ (b^c').x,
e_R3(x) = a.S^-1(x) ^ c'.x.  The second and third factors use
sum_x (-1)^{a.S^-1(x) ^ q.x} = LAT[a][q]; verify.py checks this identity and
the exhaustive t=0 correlation.)

Output (fixed format, one row per triple):
   rank  a  b  c  cprime  tau_num  tau  log2|tau|   linear-independence rank

The script writes ``masks.json`` next to itself (top list + the chosen
independent basis of 3 triples and the 7 nonzero span combinations) and prints
the same table so it can be pasted into the report.

Run:  python3.10 search_masks.py [--top N]
"""

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_tables as bt

HERE = os.path.dirname(os.path.abspath(__file__))


def gf2_rank_ints(vecs):
    basis = []
    for v in vecs:
        x = v
        while x:
            hb = x.bit_length() - 1
            found = False
            for i, b in enumerate(basis):
                if b.bit_length() - 1 == hb:
                    x ^= b
                    found = True
                    break
            if not found:
                basis.append(x)
                basis.sort(reverse=True)
                break
    return len(basis)


def encode(a, b, c):
    return (a << 16) | (b << 8) | c


def find_top(T, topn=40):
    """Return list of (tau_unit, a, b, c) sorted by |tau_unit| desc."""
    lat = T.lat.astype(np.int64)
    rows = []
    for c in range(1, 256):
        cp = bt.sigma_inv(c)
        f1 = lat[1:256, c]                                   # over b
        f3 = lat[1:256, cp]                                  # over a
        # f2[a,b] = LAT[a][b ^ cp]
        cols = np.arange(1, 256) ^ cp
        f2 = lat[1:256][:, cols]                             # (a, b)
        P = (f1[None, :] * f2) * f3[:, None]
        flat = np.abs(P).ravel()
        k = min(topn, flat.size)
        idx = np.argpartition(flat, flat.size - k)[flat.size - k:]
        for j in idx:
            a_i, b_i = divmod(int(j), 255)
            rows.append((int(P[a_i, b_i]), a_i + 1, b_i + 1, c))
    rows.sort(key=lambda r: (-abs(r[0]), r[1], r[2], r[3]))
    # de-duplicate triples (enumerated once, but be safe)
    seen = set()
    out = []
    for r in rows:
        key = (r[1], r[2], r[3])
        if key not in seen:
            seen.add(key)
            out.append(r)
    return out


def independent_basis(top, v=3):
    """Greedy: walk down the |tau| order and keep triples that keep the GF(2)
    rank growing, until v independent triples are collected."""
    chosen = []
    for row in top:
        a, b, c = row[1], row[2], row[3]
        vecs = [encode(*[x[1], x[2], x[3]]) for x in chosen] + [encode(a, b, c)]
        if gf2_rank_ints(vecs) == len(vecs):
            chosen.append(row)
            if len(chosen) == v:
                break
    return chosen


def span_masks(basis):
    """The 7 nonzero GF(2)-combinations (a,b,c) of the basis triples."""
    out = []
    v = len(basis)
    for w in range(1, 1 << v):
        a = b = c = 0
        members = []
        for i in range(v):
            if (w >> i) & 1:
                _, ai, bi, ci = basis[i]
                a ^= ai
                b ^= bi
                c ^= ci
                members.append(i)
        out.append((w, a, b, c, members))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=40)
    ap.add_argument("--basis", type=int, default=3,
                    help="dimension v of the span used for the capacity test")
    ap.add_argument("--json", default=os.path.join(HERE, "masks.json"))
    args = ap.parse_args()

    T = bt.get_tables()
    top = find_top(T, args.top)
    basis = independent_basis(top, args.basis)
    if len(basis) < args.basis:
        raise SystemExit("could not find %d linearly independent triples" % args.basis)

    print("# Reduced-scale SNOW-Axn analogue: linear masks")
    print("# tau_num = LAT[b][c]*LAT[a][b^c']*LAT[a][c'],  tau = tau_num/2^24")
    print("# sigma = rotl3, c' = sigma^-1(c) = rotr3(c)")
    print("%3s %6s %6s %6s %6s %10s %14s %12s" %
          ("rank", "a", "b", "c", "c'", "tau_num", "tau", "log2|tau|"))
    for i, (tn, a, b, c) in enumerate(top[:args.top]):
        if tn == 0:
            break
        tau = tn / 16777216.0
        l2 = (np.log2(abs(tn)) - 24.0) if tn else float("-inf")
        print("%3d  0x%02x   0x%02x   0x%02x   0x%02x  %10d  %+14.9f  %12.6f" %
              (i + 1, a, b, c, bt.sigma_inv(c), tn, tau, l2))
        if i >= 19:
            break

    print()
    print("# chosen linearly independent basis (v=%d):" % args.basis)
    for i, (tn, a, b, c) in enumerate(basis):
        tau = tn / 16777216.0
        print("  t%d: a=0x%02x b=0x%02x c=0x%02x c'=0x%02x  tau=%+0.9f  log2|tau|=%0.6f" %
              (i + 1, a, b, c, bt.sigma_inv(c), tau, np.log2(abs(tau))))
    print("# rank of the %d basis triples as GF(2)^24 vectors: %d" %
          (len(basis), gf2_rank_ints([encode(x[1], x[2], x[3]) for x in basis])))

    span = span_masks(basis)
    print()
    print("# nonzero span combinations (XOR of the basis masks)")
    print("%2s %6s %6s %6s %6s %6s %10s %14s %12s" %
          ("w", "a", "b", "c", "c'", "members", "tau_num", "tau", "log2|tau|"))
    span_rows = []
    for (w, a, b, c, members) in span:
        tn = bt.tau_unit(a, b, c, T)
        tau = tn / 16777216.0
        l2 = (np.log2(abs(tn)) - 24.0) if tn else float("-inf")
        print("%2d  0x%02x   0x%02x   0x%02x   0x%02x   %-9s %10d  %+14.9f  %12.6f" %
              (w, a, b, c, bt.sigma_inv(c),
               "+".join("t%d" % (m + 1) for m in members), tn, tau, l2))
        span_rows.append({"w": w, "a": a, "b": b, "c": c,
                          "cprime": bt.sigma_inv(c),
                          "members": members, "tau_num": tn, "tau": tau})

    pred_sei2 = sum(float(r["tau"]) ** 4 for r in span_rows)
    print()
    print("# predicted capacity SEI_2(V_3) = sum_w tau(w)^4 = %.6e  (log2 = %.6f)" %
          (pred_sei2, np.log2(pred_sei2) if pred_sei2 else float("-inf")))

    payload = {
        "sigma": "rotl3",
        "tau_formula": "LAT[b][c]*LAT[a][b^c']*LAT[a][c'] / 2^24",
        "top": [{"a": a, "b": b, "c": c, "cprime": bt.sigma_inv(c),
                 "tau_num": tn, "tau": tn / 16777216.0} for (tn, a, b, c) in top[:100]],
        "basis": [{"a": a, "b": b, "c": c, "cprime": bt.sigma_inv(c),
                   "tau_num": tn, "tau": tn / 16777216.0} for (tn, a, b, c) in basis],
        "span": span_rows,
        "predicted_SEI2": pred_sei2,
    }
    with open(args.json, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=1)
    print("# wrote", args.json)


if __name__ == "__main__":
    main()
