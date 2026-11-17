"""Profiles of the four feasible support pairs.

Inputs: the tables of verify.py.
Outputs: for each support of Table 2 the number of nonzero pairs, the number of
         levels, the maximum and the weakest level.
Paper object: Section 4.1 and the paragraph on the other three supports.
"""

import sys, io, contextlib, os, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(CODE, "mask_search", "correlation"))
with contextlib.redirect_stdout(io.StringIO()):
    from verify import LAT, ILAT, L, Linv, sigma, pullback

SUPPORTS = [(4, 5, 6, 12, 13, 14), (0, 2, 3, 8, 10, 11), (4, 6, 7, 12, 14, 15), (0, 1, 2, 8, 9, 10)]


def nullspace(rows, n):
    rows = list(rows); piv = {}; r = 0
    for col in range(n):
        idx = None
        for i in range(r, len(rows)):
            if (rows[i] >> col) & 1:
                idx = i; break
        if idx is None:
            continue
        rows[r], rows[idx] = rows[idx], rows[r]
        for i in range(len(rows)):
            if i != r and ((rows[i] >> col) & 1):
                rows[i] ^= rows[r]
        piv[col] = r; r += 1
        if r == len(rows):
            break
    out = []
    for f in [c for c in range(n) if c not in piv]:
        v = 1 << f
        for c, ri in piv.items():
            if (rows[ri] >> f) & 1:
                v |= 1 << c
        out.append(v)
    return out


def admissible_basis(mu, f):
    free = [(b, t) for b in mu for t in range(8)]
    outc = [b for b in range(16) if b not in mu]
    img = [f(1 << (8 * b + t)) for (b, t) in free]
    rows = []
    for ob in outc:
        for ot in range(8):
            row = 0
            for j, w in enumerate(img):
                if (w >> (8 * ob + ot)) & 1:
                    row |= 1 << j
            rows.append(row)
    out = []
    for v in nullspace(rows, len(free)):
        m = 0
        for j, (b, t) in enumerate(free):
            if (v >> j) & 1:
                m |= 1 << (8 * b + t)
        out.append(m)
    return out


def arr(v):
    return np.array([[x >> (8 * i) & 0xFF for i in range(16)] for x in v], dtype=np.uint8)


def expand(basis, k, sel):
    return np.bitwise_xor.reduce(sel[:, :, None] * basis[None, :, :], axis=1)


absL = np.abs(LAT.astype(np.int16))
absI = np.abs(ILAT.astype(np.int16))
ILa = np.abs(ILAT.astype(np.int32))
T12a = np.zeros((256, 256), dtype=np.uint8)
for m1 in range(256):
    pr = ILa[m1][None, :] * ILa
    pr[:, 0] = 0
    T12a[m1] = pr.argmax(axis=1).astype(np.uint8)

for mu in SUPPORTS:
    tA = admissible_basis(mu, lambda w: pullback(w, Linv))
    tC = admissible_basis(mu, lambda w: pullback(w, lambda s: Linv(sigma(L(s)))))
    Pb_ = [pullback(b, Linv) for b in tA]
    Qb_ = [pullback(c, lambda s: Linv(sigma(L(s)))) for c in tC]
    B = expand(arr(tA), 16, ((np.arange(1 << 16)[:, None] >> np.arange(16)) & 1).astype(np.uint8))
    P = expand(arr(Pb_), 16, ((np.arange(1 << 16)[:, None] >> np.arange(16)) & 1).astype(np.uint8))
    C = expand(arr(tC), 8, ((np.arange(1 << 8)[:, None] >> np.arange(8)) & 1).astype(np.uint8))
    Q = expand(arr(Qb_), 8, ((np.arange(1 << 8)[:, None] >> np.arange(8)) & 1).astype(np.uint8))
    acc = np.zeros((len(B), len(C)))
    zero = np.zeros((len(B), len(C)), dtype=bool)
    for i in mu:
        b = B[:, i][:, None].astype(np.int32)
        c = C[:, i][None, :].astype(np.int32)
        p = P[:, i][:, None].astype(np.int32)
        q = Q[:, i][None, :].astype(np.int32)
        v1 = absL[p, q].astype(np.int32)
        cb = (b ^ c).astype(np.int32)
        a = T12a[cb, c].astype(np.int32)
        v2 = absI[cb, a].astype(np.int32)
        v3 = absI[c, a].astype(np.int32)
        zero |= (v1 == 0) | (v2 == 0) | (v3 == 0)
        acc += np.log2(np.maximum(v1, 1) / 256) + np.log2(np.maximum(v2, 1) / 256) + np.log2(np.maximum(v3, 1) / 256)
    acc[zero] = -1e9
    acc[0, 0] = -1e9
    nz = acc[acc > -1e8]
    vals = np.unique(np.round(nz, 6))
    print("support %s: nonzero=%d levels=%d max=%.3f min=%.2f"
          % (mu, nz.size, vals.size, vals.max(), vals.min()))
