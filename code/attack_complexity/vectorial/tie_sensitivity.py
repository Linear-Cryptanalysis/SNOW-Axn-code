"""Tie-breaking sensitivity of the span V_26 of the vectorial attack.

The bytewise maximizer of Equation (10) is not unique on every active byte of
the 26 triples of data/w26.json.  This script rebuilds the basis with different
byte-level tie rules, verifies that the 26 triples are linearly independent in
every case, and recomputes

    SEI_2(V_26) = sum_{w != 0} tau(w)^4 ,   SEI_3(V_26) = sum_{w != 0} tau(w)^6

by exhaustive enumeration of the 2^26 elements.  The power sums give the
positivity of log2 N = log2(4 B_g ln 2) - log2 SEI_2 and of the data of the
distinguisher of Section 4.3.

Rules (B and C are kept from data/w26.json, only A changes):
  R0  the representative stored in w26.json (the smallest maximizer)
  R1  smallest maximizer of every active byte
  R2  largest maximizer
  R3  uniformly random maximizer per active byte, fixed seed

The script validates itself on the nested spans of the first v triples
(v <= 12), whose SEI_2 values are the rows of Table 5.

Inputs:  ../data/w26.json
Outputs: for every rule the GF(2) rank, log2 SEI_2, log2 SEI_3 and the changes
         of these quantities relative to R0, printed and written to
         tie_sensitivity.json in the current directory.
Paper objects: the tie sentence of Section 4.2 and the capacities quoted in
Section 4.3.
"""

import contextlib
import io
import json
import math
import os
import random
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(CODE, "mask_search", "correlation"))
with contextlib.redirect_stdout(io.StringIO()):
    from verify import LAT, ILAT, L, Linv, sigma, pullback

MU = (4, 5, 6, 12, 13, 14)
DATA = os.path.normpath(os.path.join(HERE, "..", "..", "..", "data", "w26.json"))
LO = 20
K = 26
BLK = 1 << LO
NHI = 1 << (K - LO)
SEED = 20261118

TABLE5 = {1: -219.082321, 2: -218.656208, 3: -218.327675, 4: -218.060254,
          5: -217.834737, 6: -217.639752, 7: -217.492845, 8: -217.437326,
          9: -217.383865, 10: -217.336328, 11: -217.276048, 12: -217.235331}

absL = np.abs(LAT.astype(np.int16))
absI = np.abs(ILAT.astype(np.int16))
ILa = np.abs(ILAT.astype(np.int32))
with np.errstate(divide="ignore"):
    LL = np.log2(absL.astype(np.float64) / 256.0)
    II = np.log2(absI.astype(np.float64) / 256.0)


def to_bytes(h):
    return [int(h[2 * i:2 * i + 2], 16) for i in range(16)]


def to_word(b):
    v = 0
    for i in range(16):
        v |= int(b[i]) << (8 * i)
    return int(v)


def pull_bytecols(w, fn):
    return [(pullback(w, fn) >> (8 * i)) & 0xFF for i in range(16)]


def maximizers(m1, m2):
    pr = ILa[m1][1:].astype(np.int64) * ILa[m2][1:].astype(np.int64)
    return (np.nonzero(pr == pr.max())[0] + 1)


def rank_gf2(words):
    rows = [int(x) for x in words]
    r = 0
    for col in range(383, -1, -1):
        idx = None
        for i in range(r, len(rows)):
            if (rows[i] >> col) & 1:
                idx = i
                break
        if idx is None:
            continue
        rows[r], rows[idx] = rows[idx], rows[r]
        for i in range(len(rows)):
            if i != r and ((rows[i] >> col) & 1):
                rows[i] ^= rows[r]
        r += 1
        if r == len(rows):
            break
    return r


def load_rules():
    ent = json.load(open(DATA, encoding="utf-8"))["entries"]
    assert len(ent) == K
    A0 = np.array([to_bytes(e["A"]) for e in ent], dtype=np.uint8)
    Bb = np.array([to_bytes(e["B"]) for e in ent], dtype=np.uint8)
    Cb = np.array([to_bytes(e["C"]) for e in ent], dtype=np.uint8)
    Pb = np.array([pull_bytecols(to_word(r), Linv) for r in Bb], dtype=np.uint8)
    Qb = np.array([pull_bytecols(to_word(r), lambda s: Linv(sigma(L(s))))
                   for r in Cb], dtype=np.uint8)
    rules = {"R0": A0.copy(), "R1": A0.copy(), "R2": A0.copy(), "R3": A0.copy()}
    rng = random.Random(SEED)
    tied = 0
    for k in range(K):
        for i in MU:
            cb, c = int(Cb[k][i]) ^ int(Bb[k][i]), int(Cb[k][i])
            ms = maximizers(cb, c)
            assert int(A0[k][i]) in ms
            if len(ms) > 1:
                tied += 1
            rules["R1"][k][i] = int(ms[0])
            rules["R2"][k][i] = int(ms[-1])
            rules["R3"][k][i] = int(ms[rng.randrange(len(ms))])
    return ent, A0, Bb, Cb, Pb, Qb, rules, tied


def small_span_sei2(Ab, Bb, Cb, Pb, Qb, v):
    sel = ((np.arange(1 << v)[:, None] >> np.arange(v)) & 1).astype(np.uint8)
    AA = np.bitwise_xor.reduce(sel[:, :, None] * Ab[None, :v, :], axis=1)
    BB = np.bitwise_xor.reduce(sel[:, :, None] * Bb[None, :v, :], axis=1)
    CC = np.bitwise_xor.reduce(sel[:, :, None] * Cb[None, :v, :], axis=1)
    PP = np.bitwise_xor.reduce(sel[:, :, None] * Pb[None, :v, :], axis=1)
    QQ = np.bitwise_xor.reduce(sel[:, :, None] * Qb[None, :v, :], axis=1)
    lg = np.zeros(len(sel))
    for j in MU:
        lg += LL[PP[:, j], QQ[:, j]]
        lg += II[CC[:, j] ^ BB[:, j], AA[:, j]]
        lg += II[CC[:, j], AA[:, j]]
    lg[0] = -np.inf
    return math.log2(float(np.exp2(4.0 * lg).sum()))


def main():
    t_all = time.time()
    ent, A0, Bb, Cb, Pb, Qb, rules, tied = load_rules()
    order = ["R0", "R1", "R2", "R3"]
    print("tied active bytes over the 26 triples: %d" % tied)
    print("R1 equals the stored representatives: %s" % bool(np.array_equal(rules["R1"], rules["R0"])))

    ok = True
    for v in range(1, 13):
        got = small_span_sei2(A0, Bb, Cb, Pb, Qb, v)
        ok &= abs(got - TABLE5[v]) < 1e-5
    print("validation against Table 5 for v = 1..12: %s" % ("PASS" if ok else "FAIL"))
    if not ok:
        raise SystemExit("validation failed")

    bitslo = ((np.arange(BLK, dtype=np.int64)[None, :]
               >> np.arange(LO, dtype=np.int64)[:, None]) & 1).astype(np.uint8)
    act = list(MU)

    def expand_low(bas):
        out = np.zeros((len(act), BLK), dtype=np.uint8)
        for k in range(LO):
            out ^= bas[k][:, None] * bitslo[k][None, :]
        return out

    def expand_high(bas):
        out = np.zeros((NHI, len(act)), dtype=np.uint8)
        for hi in range(1, NHI):
            v = np.zeros(len(act), dtype=np.uint8)
            for t in range(K - LO):
                if (hi >> t) & 1:
                    v ^= bas[LO + t]
            out[hi] = v
        return out

    Blo, Clo = expand_low(Bb[:, act]), expand_low(Cb[:, act])
    Plo, Qlo = expand_low(Pb[:, act]), expand_low(Qb[:, act])
    CBl = Clo ^ Blo
    Bhi, Chi = expand_high(Bb[:, act]), expand_high(Cb[:, act])
    Phi, Qhi = expand_high(Pb[:, act]), expand_high(Qb[:, act])
    Alo = {rn: expand_low(rules[rn][:, act]) for rn in order}
    Ahi = {rn: expand_high(rules[rn][:, act]) for rn in order}

    s4 = {rn: 0.0 for rn in order}
    s6 = {rn: 0.0 for rn in order}
    t0 = time.time()
    for hi in range(NHI):
        Cfull = Clo ^ Chi[hi][:, None]
        CBfull = CBl ^ (Chi[hi] ^ Bhi[hi])[:, None]
        Pfull = Plo ^ Phi[hi][:, None]
        Qfull = Qlo ^ Qhi[hi][:, None]
        lg0 = np.zeros(BLK)
        for j in range(len(act)):
            lg0 += LL[Pfull[j], Qfull[j]]
        for rn in order:
            A = Alo[rn] ^ Ahi[rn][hi][:, None]
            lg = lg0.copy()
            for j in range(len(act)):
                lg += II[CBfull[j], A[j]]
                lg += II[Cfull[j], A[j]]
            if hi == 0:
                lg[0] = -np.inf
            s4[rn] += float(np.exp2(4.0 * lg).sum())
            s6[rn] += float(np.exp2(6.0 * lg).sum())
        if hi % 16 == 15 or hi == NHI - 1:
            print("  block %2d/%d  %.1fs" % (hi + 1, NHI, time.time() - t0), flush=True)

    ref2, ref3 = -216.797094, -327.184547
    out = {"data": os.path.basename(DATA), "rules": {}}
    print("\nrule  rank  log2 SEI_2   log2 SEI_3   dSEI_2    dlog2N")
    for rn in order:
        v384 = [to_word(rules[rn][k]) | (to_word(Bb[k]) << 128) | (to_word(Cb[k]) << 256)
                for k in range(K)]
        rk = rank_gf2(v384)
        v2 = math.log2(s4[rn])
        v3 = math.log2(s6[rn])
        print("%-4s %4d  %11.6f %12.6f %+9.6f %+9.6f"
              % (rn, rk, v2, v3, v2 - ref2, -0.5 * (v2 - ref2)))
        out["rules"][rn] = dict(rank384=rk, log2_SEI2=v2, log2_SEI3=v3,
                                dSEI2=v2 - ref2, dSEI3=v3 - ref3,
                                dlog2N=-0.5 * (v2 - ref2),
                                A_bytes_diff=int((rules[rn] != A0).sum()))
    out["baseline"] = {"log2_SEI2": ref2, "log2_SEI3": ref3}
    out["seconds"] = round(time.time() - t_all, 1)
    json.dump(out, open("tie_sensitivity.json", "w", encoding="utf-8"), indent=1)
    print("total %.1fs ; wrote tie_sensitivity.json" % out["seconds"])


if __name__ == "__main__":
    main()
