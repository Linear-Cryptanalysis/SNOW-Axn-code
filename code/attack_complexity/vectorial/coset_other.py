"""Effect of the three further feasible supports on the capacity of W.

Inputs: ../data/w26.json and the tables of verify.py.
Outputs: for the strongest triple t of each further feasible support, the coset
         sum of tau^4 over V_26 + t, the new SEI_2 and the change of the data
         complexity.
Paper object: the sentence of Section 4.2 on the supports outside the profile.
"""
import sys, io, contextlib, json, math, time, os
import numpy as np

CODE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(CODE, "mask_search", "correlation"))
with contextlib.redirect_stdout(io.StringIO()):
    from verify import LAT, ILAT, L, Linv, sigma, pullback

SUPPORTS = [(0, 2, 3, 8, 10, 11), (4, 6, 7, 12, 14, 15), (0, 1, 2, 8, 9, 10)]


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


def best_triple(mu):
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
    bi, cj = np.unravel_index(np.argmax(acc), acc.shape)
    Ab = [0] * 16
    for i in mu:
        Ab[i] = int(T12a[int(C[cj][i]) ^ int(B[bi][i])][int(C[cj][i])])
    return np.array(Ab, np.uint8), B[bi], C[cj], float(acc[bi, cj])


ent = json.load(open(os.path.normpath(os.path.join(CODE, "..", "data", "w26.json")), encoding="utf-8"))["entries"]


def tb(h):
    return [int(h[2 * i:2 * i + 2], 16) for i in range(16)]


Ab_ = np.array([tb(e["A"]) for e in ent], dtype=np.uint8)
Bb_ = np.array([tb(e["B"]) for e in ent], dtype=np.uint8)
Cb_ = np.array([tb(e["C"]) for e in ent], dtype=np.uint8)
Pb__ = np.array([[pullback(int.from_bytes(bytes(r), "little"), Linv) >> (8 * i) & 0xFF for i in range(16)] for r in Bb_], dtype=np.uint8)
Qb__ = np.array([[pullback(int.from_bytes(bytes(r), "little"), lambda s: Linv(sigma(L(s)))) >> (8 * i) & 0xFF for i in range(16)] for r in Cb_], dtype=np.uint8)

SEI24 = 2.0 ** -216.797094
CH = 1 << 20
offsets = []
for mu in SUPPORTS:
    At, Bt, Ct, lg = best_triple(mu)
    Pt = np.array([pullback(int.from_bytes(bytes(Bt), "little"), Linv) >> (8 * i) & 0xFF for i in range(16)], dtype=np.uint8)
    Qt = np.array([pullback(int.from_bytes(bytes(Ct), "little"), lambda s: Linv(sigma(L(s)))) >> (8 * i) & 0xFF for i in range(16)], dtype=np.uint8)
    offsets.append((mu, At, Bt, Ct, Pt, Qt, lg))
    print("support %s: best log2|tau| = %.6f" % (mu, lg), flush=True)

sums = [0.0] * len(offsets)
t0 = time.time()
for start in range(0, 1 << 26, CH):
    idx = np.arange(start, start + CH, dtype=np.int64)
    s = ((idx[:, None] >> np.arange(26)) & 1).astype(np.uint8)
    A0 = expand(Ab_, 26, s)
    B0 = expand(Bb_, 26, s)
    C0 = expand(Cb_, 26, s)
    P0 = expand(Pb__, 26, s)
    Q0 = expand(Qb__, 26, s)
    for k, (mu, At, Bt, Ct, Pt, Qt, lg) in enumerate(offsets):
        A2 = A0 ^ At
        B2 = B0 ^ Bt
        C2 = C0 ^ Ct
        P2 = P0 ^ Pt
        Q2 = Q0 ^ Qt
        lgt = np.zeros(CH)
        zr = np.zeros(CH, dtype=bool)
        for l in range(16):
            v1 = absL[P2[:, l], Q2[:, l]].astype(np.float64)
            v2 = absI[C2[:, l] ^ B2[:, l], A2[:, l]].astype(np.float64)
            v3 = absI[C2[:, l], A2[:, l]].astype(np.float64)
            zr |= (v1 == 0) | (v2 == 0) | (v3 == 0)
            lgt += np.log2(np.maximum(v1, 1) / 256) + np.log2(np.maximum(v2, 1) / 256) + np.log2(np.maximum(v3, 1) / 256)
        sums[k] += float(np.exp2(4.0 * lgt[~zr]).sum())
    if (start // CH) % 8 == 0:
        print("chunk %d/64  %.0fs" % (start // CH, time.time() - t0), flush=True)

for k, (mu, At, Bt, Ct, Pt, Qt, lg) in enumerate(offsets):
    new = SEI24 + sums[k]
    print("support %s: coset sum4 = 2^%.6f, SEI_2 becomes 2^%.6f, data change %.4f bits"
          % (mu, math.log2(sums[k]), math.log2(new), -math.log2(new / SEI24)))
print("total time %.0fs" % (time.time() - t0))

