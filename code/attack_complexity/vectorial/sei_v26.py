"""Power sums of the span of the 26 strongest linear mask triples.

Inputs:  ../data/w24.json (22 highest triples plus two of level 8) and
         ../data/w26.json (the same 22 plus all four triples of level 8).
Outputs: SEI_2(V_24), SEI_2(V_25), SEI_2(V_26) and SEI_3(V_24), SEI_3(V_26).

V_26 is split into the cosets V_24, V_24 + t1, V_24 + t2 and V_24 + t1 + t2,
where t1 and t2 are the two triples of w26.json that are not in w24.json.
Paper object: Table 5 rows v = 25 and v = 26, the data complexity of
Section 4.2 and the distinguisher capacity of Section 4.3.
"""

import sys, io, contextlib, json, math, os, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(CODE, "mask_search", "correlation"))
with contextlib.redirect_stdout(io.StringIO()):
    from verify import LAT, ILAT, L, Linv, sigma, pullback

MU = (4, 5, 6, 12, 13, 14)
DATA24 = os.path.normpath(os.path.join(CODE, "..", "data", "w24.json"))
DATA26 = os.path.normpath(os.path.join(CODE, "..", "data", "w26.json"))


def to_bytes(h):
    return [int(h[2 * i:2 * i + 2], 16) for i in range(16)]


d24 = json.load(open(DATA24, encoding="utf-8"))["entries"]
d26 = json.load(open(DATA26, encoding="utf-8"))["entries"]
base_keys = {e["A"] + e["B"] + e["C"] for e in d24}
extra = [e for e in d26 if e["A"] + e["B"] + e["C"] not in base_keys]
assert len(extra) == 2

Ab_ = np.array([to_bytes(e["A"]) for e in d24], dtype=np.uint8)
Bb_ = np.array([to_bytes(e["B"]) for e in d24], dtype=np.uint8)
Cb_ = np.array([to_bytes(e["C"]) for e in d24], dtype=np.uint8)
Pb__ = np.array([[pullback(int.from_bytes(bytes(r), "little"), Linv) >> (8 * i) & 0xFF for i in range(16)] for r in Bb_], dtype=np.uint8)
Qb__ = np.array([[pullback(int.from_bytes(bytes(r), "little"), lambda s: Linv(sigma(L(s)))) >> (8 * i) & 0xFF for i in range(16)] for r in Cb_], dtype=np.uint8)


def offset(e):
    At = np.array(to_bytes(e["A"]), dtype=np.uint8)
    Bt = np.array(to_bytes(e["B"]), dtype=np.uint8)
    Ct = np.array(to_bytes(e["C"]), dtype=np.uint8)
    Pt = np.array([pullback(int.from_bytes(bytes(Bt), "little"), Linv) >> (8 * i) & 0xFF for i in range(16)], dtype=np.uint8)
    Qt = np.array([pullback(int.from_bytes(bytes(Ct), "little"), lambda s: Linv(sigma(L(s)))) >> (8 * i) & 0xFF for i in range(16)], dtype=np.uint8)
    return At, Bt, Ct, Pt, Qt


o1, o2 = offset(extra[0]), offset(extra[1])
comb = tuple(a ^ b for a, b in zip(o1, o2))
variants = [None, o1, o2, comb]


def expand(basis, k, sel):
    return np.bitwise_xor.reduce(sel[:, :, None] * basis[None, :, :], axis=1)


absL = np.abs(LAT.astype(np.int16))
absI = np.abs(ILAT.astype(np.int16))

t0 = time.time()
CH = 1 << 20
sums4 = [0.0] * 4
sums6 = [0.0] * 4
for start in range(0, 1 << 24, CH):
    idx = np.arange(start, start + CH, dtype=np.int64)
    s = ((idx[:, None] >> np.arange(24)) & 1).astype(np.uint8)
    A = expand(Ab_, 24, s)
    Bc = expand(Bb_, 24, s)
    Cc = expand(Cb_, 24, s)
    Pp = expand(Pb__, 24, s)
    Qq = expand(Qb__, 24, s)
    for vi, v in enumerate(variants):
        if v is None:
            A2, B2, C2, P2, Q2 = A, Bc, Cc, Pp, Qq
        else:
            A2, B2, C2, P2, Q2 = A ^ v[0], Bc ^ v[1], Cc ^ v[2], Pp ^ v[3], Qq ^ v[4]
        lg = np.zeros(CH)
        zr = np.zeros(CH, dtype=bool)
        for l in range(16):
            v1 = absL[P2[:, l], Q2[:, l]].astype(np.float64)
            v2 = absI[C2[:, l] ^ B2[:, l], A2[:, l]].astype(np.float64)
            v3 = absI[C2[:, l], A2[:, l]].astype(np.float64)
            zr |= (v1 == 0) | (v2 == 0) | (v3 == 0)
            lg += np.log2(np.maximum(v1, 1) / 256) + np.log2(np.maximum(v2, 1) / 256) + np.log2(np.maximum(v3, 1) / 256)
        if start == 0 and vi == 0:
            zr[0] = True
        sums4[vi] += float(np.exp2(4.0 * lg[~zr]).sum())
        sums6[vi] += float(np.exp2(6.0 * lg[~zr]).sum())

s4_24 = sums4[0]
s4_25 = sums4[0] + sums4[1]
s4_26 = sums4[0] + sums4[1] + sums4[2] + sums4[3]
s6_24 = sums6[0]
s6_26 = sums6[0] + sums6[1] + sums6[2] + sums6[3]
print("SEI_2(V_24) = 2^%.6f" % math.log2(s4_24))
print("SEI_2(V_25) = 2^%.6f" % math.log2(s4_25))
print("SEI_2(V_26) = 2^%.6f" % math.log2(s4_26))
print("SEI_3(V_24) = 2^%.6f" % math.log2(s6_24))
print("SEI_3(V_26) = 2^%.6f" % math.log2(s6_26))
print("log2 N at v=26 = %.4f" % (math.log2(4 * 214 * math.log(2)) - math.log2(s4_26)))
print("%.1fs" % (time.time() - t0))
