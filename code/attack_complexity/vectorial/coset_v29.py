"""Capacity of V_29: add the 29th strongest triple to V_28.

Inputs: ../data/w26.json and the tables of verify.py.
Outputs: the coset sum of tau^4 over V_28 + t, the new SEI_2 and the change of
         the data complexity.
Paper object: the joint-scan sentence of Section 4.2.
"""
import sys, io, contextlib, json, math, time, os
import numpy as np

CODE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(CODE, "mask_search", "correlation"))
with contextlib.redirect_stdout(io.StringIO()):
    from verify import LAT, ILAT, L, Linv, sigma, pullback

MU = (4, 5, 6, 12, 13, 14)


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

tA = admissible_basis(MU, lambda w: pullback(w, Linv))
tC = admissible_basis(MU, lambda w: pullback(w, lambda s: Linv(sigma(L(s)))))
Pb_ = [pullback(b, Linv) for b in tA]
Qb_ = [pullback(c, lambda s: Linv(sigma(L(s)))) for c in tC]
B = expand(arr(tA), 16, ((np.arange(1 << 16)[:, None] >> np.arange(16)) & 1).astype(np.uint8))
P = expand(arr(Pb_), 16, ((np.arange(1 << 16)[:, None] >> np.arange(16)) & 1).astype(np.uint8))
C = expand(arr(tC), 8, ((np.arange(1 << 8)[:, None] >> np.arange(8)) & 1).astype(np.uint8))
Q = expand(arr(Qb_), 8, ((np.arange(1 << 8)[:, None] >> np.arange(8)) & 1).astype(np.uint8))
acc = np.zeros((len(B), len(C)))
zero = np.zeros((len(B), len(C)), dtype=bool)
for i in MU:
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
order = np.argsort(-acc.ravel())
top = []
for k in range(32):
    bi, cj = np.unravel_index(order[k], acc.shape)
    Ab = [0] * 16
    for i in MU:
        Ab[i] = int(T12a[int(C[cj][i]) ^ int(B[bi][i])][int(C[cj][i])])
    top.append((np.array(Ab, np.uint8), B[bi], C[cj], float(acc[bi, cj])))
print("29th triple: log2|tau| = %.6f" % top[28][3], flush=True)
for k in (26, 27, 28):
    A, Bt, Ct, lg = top[k]
    print("  #%d log2|tau| = %.6f A=%s B=%s C=%s" % (k+1, lg,
          "".join("%02x" % x for x in A), "".join("%02x" % x for x in Bt), "".join("%02x" % x for x in Ct)), flush=True)

ent = json.load(open(os.path.normpath(os.path.join(CODE, "..", "data", "w26.json")), encoding="utf-8"))["entries"]
def tb(h):
    return [int(h[2 * i:2 * i + 2], 16) for i in range(16)]

A_basis = [tb(e["A"]) for e in ent]
B_basis = [tb(e["B"]) for e in ent]
C_basis = [tb(e["C"]) for e in ent]
for k in (26, 27):
    A_basis.append(list(top[k][0])); B_basis.append(list(top[k][1])); C_basis.append(list(top[k][2]))
A_basis = np.array(A_basis, dtype=np.uint8)
B_basis = np.array(B_basis, dtype=np.uint8)
C_basis = np.array(C_basis, dtype=np.uint8)


def pull(hB):
    Xi = pullback(int.from_bytes(bytes(hB), "little"), Linv)
    return np.array([(Xi >> (8 * i)) & 0xFF for i in range(16)], dtype=np.uint8)


def pullq(hC):
    Xi = pullback(int.from_bytes(bytes(hC), "little"), lambda s: Linv(sigma(L(s))))
    return np.array([(Xi >> (8 * i)) & 0xFF for i in range(16)], dtype=np.uint8)


P_basis = np.array([pull(b) for b in B_basis], dtype=np.uint8)
Q_basis = np.array([pullq(c) for c in C_basis], dtype=np.uint8)
At = np.array(top[28][0], dtype=np.uint8)
Bt = np.array(top[28][1], dtype=np.uint8)
Ct = np.array(top[28][2], dtype=np.uint8)
Pt = pull(top[28][1])
Qt = pullq(top[28][2])
SEI28 = 2.0 ** -216.741826
K = 28
CH = 1 << 20
sum4 = 0.0
t0 = time.time()
for start in range(0, 1 << K, CH):
    idx = np.arange(start, start + CH, dtype=np.int64)
    s = ((idx[:, None] >> np.arange(K)) & 1).astype(np.uint8)
    A2 = expand(A_basis, K, s) ^ At
    B2 = expand(B_basis, K, s) ^ Bt
    C2 = expand(C_basis, K, s) ^ Ct
    P2 = expand(P_basis, K, s) ^ Pt
    Q2 = expand(Q_basis, K, s) ^ Qt
    lgt = np.zeros(CH)
    zr = np.zeros(CH, dtype=bool)
    for l in range(16):
        v1 = absL[P2[:, l], Q2[:, l]].astype(np.float64)
        v2 = absI[C2[:, l] ^ B2[:, l], A2[:, l]].astype(np.float64)
        v3 = absI[C2[:, l], A2[:, l]].astype(np.float64)
        zr |= (v1 == 0) | (v2 == 0) | (v3 == 0)
        lgt += np.log2(np.maximum(v1, 1) / 256) + np.log2(np.maximum(v2, 1) / 256) + np.log2(np.maximum(v3, 1) / 256)
    sum4 += float(np.exp2(4.0 * lgt[~zr]).sum())
    if (start // CH) % 32 == 0:
        print("chunk %d/%d  %.0fs" % (start // CH, (1 << K) // CH, time.time() - t0), flush=True)
new = SEI28 + sum4
print("coset sum4 = 2^%.6f" % math.log2(sum4))
print("SEI_2(V_29) = 2^%.6f, gain %.4f bits over V_28" % (math.log2(new), math.log2(new / SEI28)))
print("total %.0fs" % (time.time() - t0))

