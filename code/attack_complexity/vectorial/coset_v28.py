"""Capacity of V_28: add the 28th strongest triple to V_27.

Inputs: ../data/w26.json and the tables of verify.py.
Outputs: the coset sum of tau^4 over V_27 + t, the new SEI_2 and the change of
         the data complexity.
Paper object: the joint-scan sentence of Section 4.2.
"""
import sys, io, contextlib, json, math, time, os
import numpy as np

CODE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(CODE, "mask_search", "correlation"))
with contextlib.redirect_stdout(io.StringIO()):
    from verify import LAT, ILAT, L, Linv, sigma, pullback

# 27th and 28th strongest triples of the optimal support pair (level 9)
T27 = ("00000000ee1ca50000000000fa1cf500", "00000000a500f200000000003300cf00", "000000007bd6ce00000000007bd6ce00")
T28 = ("00000000fa1cf50000000000ee1ca500", "000000003300cf0000000000a500f200", "000000007bd6ce00000000007bd6ce00")

absL = np.abs(LAT.astype(np.int16))
absI = np.abs(ILAT.astype(np.int16))


def tb(h):
    return [int(h[2 * i:2 * i + 2], 16) for i in range(16)]


def expand(basis, k, sel):
    return np.bitwise_xor.reduce(sel[:, :, None] * basis[None, :, :], axis=1)


ent = json.load(open(os.path.normpath(os.path.join(CODE, "..", "data", "w26.json")), encoding="utf-8"))["entries"]
A_basis = [tb(e["A"]) for e in ent]
B_basis = [tb(e["B"]) for e in ent]
C_basis = [tb(e["C"]) for e in ent]
for h in (T27,):
    A_basis.append(tb(h[0])); B_basis.append(tb(h[1])); C_basis.append(tb(h[2]))
A_basis = np.array(A_basis, dtype=np.uint8)
B_basis = np.array(B_basis, dtype=np.uint8)
C_basis = np.array(C_basis, dtype=np.uint8)


def pull(hB):
    Bi = int.from_bytes(bytes(hB), "little")
    Xi = pullback(Bi, Linv)
    return np.array([(Xi >> (8 * i)) & 0xFF for i in range(16)], dtype=np.uint8)


def pullq(hC):
    Ci = int.from_bytes(bytes(hC), "little")
    Xi = pullback(Ci, lambda s: Linv(sigma(L(s))))
    return np.array([(Xi >> (8 * i)) & 0xFF for i in range(16)], dtype=np.uint8)


P_basis = np.array([pull(b) for b in B_basis], dtype=np.uint8)
Q_basis = np.array([pullq(c) for c in C_basis], dtype=np.uint8)
At = np.array(tb(T28[0]), dtype=np.uint8)
Pt = pull(tb(T28[1]))
Qt = pullq(tb(T28[2]))
SEI27 = 2.0 ** -216.770773
K = 27
CH = 1 << 20
sum4 = 0.0
t0 = time.time()
for start in range(0, 1 << K, CH):
    idx = np.arange(start, start + CH, dtype=np.int64)
    s = ((idx[:, None] >> np.arange(K)) & 1).astype(np.uint8)
    A2 = expand(A_basis, K, s) ^ At
    B2 = expand(B_basis, K, s) ^ np.array(tb(T28[1]), dtype=np.uint8)
    C2 = expand(C_basis, K, s) ^ np.array(tb(T28[2]), dtype=np.uint8)
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
    if (start // CH) % 16 == 0:
        print("chunk %d/%d  %.0fs" % (start // CH, (1 << K) // CH, time.time() - t0), flush=True)
new = SEI27 + sum4
print("coset sum4 = 2^%.6f" % math.log2(sum4))
print("SEI_2(V_28) = 2^%.6f, gain %.4f bits over V_27" % (math.log2(new), math.log2(new / SEI27)))
print("total %.0fs" % (time.time() - t0))

