#!/usr/bin/env python3.10
"""Exact (zero Monte-Carlo error) folded correlations of the analogue.

Why an exact computation is possible
------------------------------------
Write F(s) for the register form of e at clock t (it depends only on
s = (R1,R2,R3) at that clock).  The FSM state is uniform given the LFSR, so
with s_t uniform,

    C2(d) = E[ (-1)^{F(s_t) ^ F(s_{t+d})} ]
          = 2^-24 sum_s A0(s) A_d(s),
    A_0(s) = (-1)^{F(s)},
    A_{k+1}(s) = E_h[ A_k(Phi(s,h)) ],        h = H(t+k+1) uniform,

because a one-clock FSM step Phi(s,h) needs only the single LFSR byte h, and
the H bytes of the measured windows are iid uniform: verify.py shows that
(H(t+1),...,H(t+3)) span 8d independent linear forms of the 64-bit LFSR state
for every window used here, hence they are uniform on 2^{8d}.

The step is cheap because for fixed s only R1' = sigma(R2^R3^h) varies with h:

    A_{k+1}(R1,R2,R3) = (1/256) sum_w A_k(w, S(R1), S(R2)).

Consequently A_1(s) = T0(S(R1),S(R2))/256 with T0(r2,r3)=sum_{r1}A0(r1,r2,r3);
A_2(s) = mean_w B1(w, S(R1)) which depends on R1 only; A_3 = mean_u A_2(u,...)
is constant, and the constant equals m0 = 2^-24 sum_s A0(s).  Hence

    C2(1) = 2^-24 sum_{r1,r2} B1(r1,r2) W(r1,r2),
    C2(2) = 2^-24 sum_{r1} A_2(r1) V(r1),
    C2(d) = tau^2 exactly for d = 3 (and d = 4),

where W(r1,r2)=sum_{r3}A0, V(r1)=sum_{r2,r3}A0.  (The identity for d >= 3 uses
that the d LFSR bytes H(t+1..t+d) are iid uniform, which verify.py confirms
for every H-window of length at most 4.  For d >= 5 the toy LFSR has
linearly dependent H bytes, so the recursion does not apply and those
offsets are not measured.)  For d = 1 there is also a closed form, printed as
a cross-check:

    C2(1) = B1 * B2 * B3 * T3,
    B1 = LAT[a^b][b^c^c']/2^8,  B2 = LAT[b][c]/2^8,  B3 = LAT[a][c']/2^8,
    T3 = 2^-8 sum_y (-1)^{(a^b^c').y ^ a.S^-1(y) ^ c'.S(y)}.

--check additionally validates the collapse against brute-force averages over
the LFSR bytes h for random register states.

Run:  python3.10 exact_folded.py [--json out.json] [--check]
"""

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_tables as bt
import monte_carlo as mc

HERE = os.path.dirname(os.path.abspath(__file__))


def state_arrays(T, a, b, c):
    """A0(s)=(-1)^{F(s)} on all 2^24 register states plus the collapsed sums."""
    r1 = np.arange(256, dtype=np.uint8).reshape(256, 1, 1)
    r2 = np.arange(256, dtype=np.uint8).reshape(1, 256, 1)
    r3 = np.arange(256, dtype=np.uint8).reshape(1, 1, 256)
    A = T.inv_sbox[r2] ^ T.inv_sbox[r3]          # (1,256,256)
    B = r1 ^ r2                                  # (256,256,1)
    C = T.rotl[r2 ^ r3] ^ T.sbox[r1]             # (256,256,256)
    e = T.parity[(a & A) ^ (b & B) ^ (c & C)]
    sgn = (1 - 2 * e).astype(np.int8)
    m0 = float(sgn.mean())
    W = sgn.sum(axis=2, dtype=np.int32)          # [r1,r2]
    V = sgn.sum(axis=(1, 2), dtype=np.int32)     # [r1]
    T0 = sgn.sum(axis=0, dtype=np.int32)         # [r2,r3]
    sbox = T.sbox.astype(np.int64)
    B1 = T0[sbox[:, None], sbox[None, :]] / 256.0
    A2 = B1.mean(axis=0)[sbox]                   # A_2(s), function of r1 only
    return {"sgn": sgn, "W": W, "V": V, "T0": T0, "B1": B1, "A2": A2,
            "m0": m0, "D3": float(A2.mean())}


def closed_form_C2_1(T, a, b, c):
    cp = bt.sigma_inv(c)
    lat = T.lat.astype(np.int64)
    y = np.arange(256, dtype=np.int64)
    sgn3 = (1 - 2 * (T.parity[((a ^ b ^ cp) & y)] ^
                     T.parity[(a & T.inv_sbox.astype(np.int64))] ^
                     T.parity[(cp & T.sbox.astype(np.int64))]).astype(np.int64))
    T3 = float(np.mean(sgn3))
    B1 = lat[(a ^ b) & 255][(b ^ c ^ cp) & 255] / 256.0
    B2 = lat[b][c] / 256.0
    B3 = lat[a][cp] / 256.0
    return B1 * B2 * B3 * T3, {"B1": B1, "B2": B2, "B3": B3, "T3": T3}


def exact_values(T, a, b, c):
    arr = state_arrays(T, a, b, c)
    C2_1 = float((arr["B1"] * arr["W"]).sum() / float(1 << 24))
    C2_2 = float((arr["A2"] * arr["V"]).sum() / float(1 << 24))
    C2_3 = arr["D3"] * arr["m0"]
    cf, factors = closed_form_C2_1(T, a, b, c)
    return {"m0": arr["m0"], "C2_1": C2_1, "C2_2": C2_2, "C2_3": C2_3,
            "C2_1_closed": cf, "factors": factors, "arrays": arr}


def brute_force_checks(T, a, b, c, arr, rng, nstates=12):
    """Validate the collapsed arrays against direct enumeration over LFSR bytes.

    A_1(s) = mean_h A0(Phi(s,h))              (256 evaluations)
    A_2(s) = mean_{h1,h2} A0(Phi(Phi(s,h1),h2))   (65536 evaluations)
    """
    sgn = arr["sgn"]
    sbox = T.sbox.astype(np.int64)
    rotl = T.rotl.astype(np.int64)
    maxdev1 = 0.0
    maxdev2 = 0.0
    for _ in range(nstates):
        r1, r2, r3 = (int(v) for v in rng.integers(0, 256, size=3))
        # A_1
        tot = 0.0
        for h in range(256):
            tot += sgn[rotl[r2 ^ r3 ^ h], sbox[r1], sbox[r2]]
        a1 = tot / 256.0
        exp1 = arr["T0"][sbox[r1], sbox[r2]] / 256.0
        maxdev1 = max(maxdev1, abs(a1 - exp1))
        # A_2
        tot = 0.0
        for h1 in range(256):
            u1 = rotl[r2 ^ r3 ^ h1]
            v1 = sbox[r1]
            w1 = sbox[r2]
            for h2 in range(256):
                tot += sgn[rotl[v1 ^ w1 ^ h2], sbox[u1], sbox[v1]]
        a2 = tot / 65536.0
        exp2 = arr["A2"][r1]
        maxdev2 = max(maxdev2, abs(a2 - exp2))
    return maxdev1, maxdev2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default=None)
    ap.add_argument("--check", action="store_true",
                    help="validate the collapsed arrays by brute force")
    args = ap.parse_args()
    T = bt.get_tables()
    basis, _ = mc.load_basis()
    combos = mc.combo_masks(basis)
    print("# Exact folded correlations of the reduced-scale analogue")
    print("# (zero Monte-Carlo error; see module docstring for the derivation)")
    print("# w    a      b      c      c'   | tau_pred     tau_exact    "
          "C2(1)        C2(2)        C2(d>=3)     C2(1)/tau^2   bits")
    out = {"combos": [], "basis": [list(x) for x in basis]}
    sei_exact = 0.0
    sei_pred = 0.0
    rng = np.random.default_rng(20261121)
    for (w, a, b, c) in combos:
        res = exact_values(T, a, b, c)
        tau = bt.tau_unit(a, b, c, T) / float(1 << 24)
        c2 = res["C2_1"]
        ratio = c2 / (tau * tau) if tau else float("nan")
        bits = (np.log2(abs(c2)) - 2 * np.log2(abs(tau))) if (tau and c2) else float("nan")
        print("%d 0x%02x 0x%02x 0x%02x 0x%02x | %+0.6e %+0.6e %+0.6e "
              "%+0.6e %+0.6e %12.4f %+8.3f"
              % (w, a, b, c, bt.sigma_inv(c), tau, res["m0"], res["C2_1"],
                 res["C2_2"], res["C2_3"], ratio, bits))
        print("        closed-form C2(1) = %+0.6e  (factors B1=%+.6f B2=%+.6f "
              "B3=%+.6f T3=%+.6f, diff=%.2e)"
              % (res["C2_1_closed"], res["factors"]["B1"], res["factors"]["B2"],
                 res["factors"]["B3"], res["factors"]["T3"],
                 res["C2_1"] - res["C2_1_closed"]))
        if args.check:
            d1, d2 = brute_force_checks(T, a, b, c, res["arrays"], rng)
            print("        brute-force check: max|A_1 mean_h - closed| = %.2e, "
                  "max|A_2 mean_hh - closed| = %.2e" % (d1, d2))
        sei_exact += res["C2_1"] ** 2
        sei_pred += tau ** 4
        out["combos"].append({
            "w": w, "a": a, "b": b, "c": c, "cprime": bt.sigma_inv(c),
            "tau_pred": tau, "tau_exact": res["m0"],
            "C2_1": res["C2_1"], "C2_2": res["C2_2"], "C2_3": res["C2_3"],
            "C2_1_closed_formula": res["C2_1_closed"],
            "factors": res["factors"], "ratio_C2_to_tau2": ratio,
            "bits_C2_vs_tau2": bits,
        })
    print()
    print("# exact capacity SEI_2 = sum_w C2(w)^2 (d=1) = %.9e  log2 = %.6f"
          % (sei_exact, np.log2(sei_exact)))
    print("# predicted capacity SEI_2 = sum_w tau(w)^4 = %.9e  log2 = %.6f"
          % (sei_pred, np.log2(sei_pred)))
    print("# exact / predicted = %.6e  (log2 difference = %+.3f bits)"
          % (sei_exact / sei_pred, np.log2(sei_exact / sei_pred)))
    out["capacity"] = {"SEI2_exact": sei_exact, "SEI2_pred": sei_pred,
                       "log2_exact": float(np.log2(sei_exact)),
                       "log2_pred": float(np.log2(sei_pred)),
                       "log2_ratio": float(np.log2(sei_exact / sei_pred))}
    # capacity if the folded offset is d = 3 (disjoint windows; C2 = tau^2)
    sei_d3 = sum((bt.tau_unit(a, b, c, T) / float(1 << 24)) ** 4
                 for (w, a, b, c) in combos)
    print("# capacity at folded offset d=3 (C2(d)=tau^2 exactly): SEI_2 = "
          "%.9e  log2 = %.6f  (identical to the prediction)"
          % (sei_d3, np.log2(sei_d3)))
    out["capacity"]["SEI2_d3"] = sei_d3
    if args.json:
        with open(args.json, "w", encoding="utf-8") as f:
            json.dump(out, f, indent=1)
        print("# wrote", args.json)


if __name__ == "__main__":
    main()
