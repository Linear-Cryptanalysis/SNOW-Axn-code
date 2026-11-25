"""Number of linear approximations in the correlation profile of Section 4.1.

The profile of Section 4.1 records one representative for each pair (B,C), so
its counts are numbers of pairs.  The necessary condition of Theorem 4 of
Zhang et al. counts linear approximations: a pair with n_i maximizers on each
of its six active bytes corresponds to prod_i n_i triples (A,B,C) of equal
correlation.  This script sums that product over the nonzero pairs of the four
feasible activity-pattern supports.

Inputs:  none (the supports are the four rows of Table 2).
Outputs: for every support the number of nonzero pairs, the number of
         approximations and log2 of that number, the four-support totals, and
         the test of the necessary condition sum_i L_i - v > 2^v at v = 26.

Paper objects: the counts of Section 4.1 and the sentence of Section 4.2 on
the necessary condition of Theorem 4.
"""

import contextlib
import io
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
with contextlib.redirect_stdout(io.StringIO()):
    from verify import LAT, ILAT, L, Linv, sigma, pullback

MU = (4, 5, 6, 12, 13, 14)
SUPPORTS = [(4, 5, 6, 12, 13, 14),
            (0, 2, 3, 8, 10, 11),
            (4, 6, 7, 12, 14, 15),
            (0, 1, 2, 8, 9, 10)]
V = 26


def nullspace(rows, n):
    rows = list(rows)
    piv = {}
    r = 0
    for col in range(n):
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
        piv[col] = r
        r += 1
        if r == len(rows):
            break
    free = [c for c in range(n) if c not in piv]
    out = []
    for f in free:
        v = 1 << f
        for c, ri in piv.items():
            if (rows[ri] >> f) & 1:
                v |= 1 << c
        out.append(v)
    return out, piv


def admissible_basis(mu, f):
    """Basis of {m : supp(m) <= mu, f(m) vanishes off mu}, as 128-bit masks."""
    free = [(b, t) for b in mu for t in range(8)]
    outc = [b for b in range(16) if b not in mu]
    n = len(free)
    img = [f(1 << (8 * b + t)) for (b, t) in free]
    rows = []
    for ob in outc:
        for ot in range(8):
            row = 0
            for j, w in enumerate(img):
                if (w >> (8 * ob + ot)) & 1:
                    row |= 1 << j
            rows.append(row)
    ker, _ = nullspace(rows, n)
    basis = []
    for v in ker:
        mask = 0
        for j, (b, t) in enumerate(free):
            if (v >> j) & 1:
                mask |= 1 << (8 * b + t)
        basis.append(mask)
    return basis


def expand(basis, k):
    sel = ((np.arange(1 << k)[:, None] >> np.arange(k)) & 1).astype(np.uint8)
    return np.bitwise_xor.reduce(sel[:, :, None] * basis[None, :, :], axis=1)


LATi = LAT.astype(np.int16)
ILATi = ILAT.astype(np.int16)
IL = np.abs(ILATi.astype(np.int32))
# T12v: max over u != 0 of |ILAT[m1][u]| |ILAT[m2][u]|; T12c: number of maximizers
T12v = np.zeros((256, 256), dtype=np.int32)
T12c = np.zeros((256, 256), dtype=np.int16)
for m1 in range(256):
    prod = IL[m1][None, :] * IL
    prod[:, 0] = 0
    mx = prod.max(axis=1)
    T12v[m1] = mx
    T12c[m1] = (prod == mx[:, None]).sum(axis=1).astype(np.int16)

t0 = time.time()
results = []
for mu in SUPPORTS:
    tA = admissible_basis(mu, lambda w: pullback(w, Linv))
    tC = admissible_basis(mu, lambda w: pullback(w, lambda s: Linv(sigma(L(s)))))
    basisB = np.array([[(v >> (8 * i)) & 0xFF for i in range(16)] for v in tA], dtype=np.uint8)
    basisP = np.array([[(pullback(v, Linv) >> (8 * i)) & 0xFF for i in range(16)]
                       for v in tA], dtype=np.uint8)
    basisC = np.array([[(v >> (8 * i)) & 0xFF for i in range(16)] for v in tC], dtype=np.uint8)
    basisQ = np.array([[(pullback(v, lambda s: Linv(sigma(L(s)))) >> (8 * i)) & 0xFF
                        for i in range(16)] for v in tC], dtype=np.uint8)
    B = expand(basisB, len(tA))
    P = expand(basisP, len(tA))
    C = expand(basisC, len(tC))
    Q = expand(basisQ, len(tC))
    zero = np.zeros((len(B), len(C)), dtype=bool)
    prodn = np.ones((len(B), len(C)), dtype=np.int64)
    for i in mu:
        b = B[:, i][:, None].astype(np.int32)
        c = C[:, i][None, :].astype(np.int32)
        p = P[:, i][:, None].astype(np.int32)
        q = Q[:, i][None, :].astype(np.int32)
        v1 = np.abs(LATi[p, q]).astype(np.int32)
        cb = (b ^ c).astype(np.int32)
        v2 = T12v[cb, c]
        zero |= (v1 == 0) | (v2 == 0)
        prodn *= np.maximum(T12c[cb, c], 1)
    zero[0, 0] = True
    nz = ~zero
    npairs = int(nz.sum())
    triples = int(prodn[nz].sum())
    results.append((mu, npairs, triples))
    print("support %s: pairs %d, approximations %d (log2 %.6f)"
          % (mu, npairs, triples, np.log2(triples)))

print()
print("pairs, optimal support: %d  (paper: 11,093,111)" % results[0][1])
print("approximations, optimal support: %d = 2^%.6f"
      % (results[0][2], np.log2(results[0][2])))
print("pairs, four supports: %d" % sum(r[1] for r in results))
print("approximations, four supports: %d = 2^%.6f"
      % (sum(r[2] for r in results), np.log2(sum(r[2] for r in results))))
print("necessary condition of Theorem 4 at v = %d: %d > %d, holds: %s"
      % (V, results[0][2] - V, 2 ** V, results[0][2] - V > 2 ** V))
print("%.1fs" % (time.time() - t0))
