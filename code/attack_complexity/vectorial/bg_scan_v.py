"""Joint scan over v and B_g for the vectorized attack.

Inputs: the capacities of V_24, V_25, V_26 (see cost.py and sei_v26.py) and of
        V_27 (see coset_v27.py).
Outputs: for each v the optimal B_g, the data, the time and the memory of the
         full attack, with the two later passes of Section 3.4 included.
Paper object: the joint-scan sentence of Section 4.2.
"""

import math
import os

PS = 0.999992
L = 512
TAU = 2401 / 2 ** 66
import json
HERE = os.path.dirname(os.path.abspath(__file__))
sei = (json.load(open(os.path.join(HERE, "sei_nested.json"))) +
       json.load(open(os.path.join(HERE, "sei_nested_24.json"))) +
       [-216.820106, -216.797094, -216.770773, -216.741826, -216.726807])
SEI = {v: s for v, s in enumerate(sei, start=1)}


def qinv(p):
    lo, hi = 0.0, 60.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if 0.5 * math.erfc(mid / math.sqrt(2)) > p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


q1 = qinv(1 - PS)


def Meq(corr, B):
    return (1 / corr ** 2) * (q1 + qinv(2.0 ** (-B - 1))) ** 2


def dec_cost(B, M):
    if B <= math.log2(M):
        return B * 2.0 ** B + M
    return 2.0 ** B * (1 + math.ceil(math.log2(M)))


_CACHE = {}

def later_passes(rem):
    if rem in _CACHE:
        return _CACHE[rem]
    best = None
    for B2 in range(rem + 1):
        B3 = rem - B2
        M2 = Meq(TAU * TAU, B2)
        N2 = 2 * M2
        T2 = N2 + dec_cost(B2, M2)
        if B3 > 0:
            M3 = Meq(TAU, B3)
            N3 = M3
            T3 = N3 + dec_cost(B3, M3)
        else:
            N3, T3 = 0.0, 0.0
        data = max(N2, N3)
        time = T2 + T3
        if best is None or max(time, data) < best[0]:
            best = (max(time, data), B2, B3, data, time)
    _CACHE[rem] = best
    return best


for v, s in SEI.items():
    S = 2.0 ** s
    rows = []
    for Bg in range(1, L + 1):
        N1 = 4 * Bg * math.log(2) / S
        T1 = 2 * N1 + Bg * 2.0 ** Bg
        _, B2, B3, dl, tl = later_passes(L - Bg)
        data = max(N1, dl)
        time = T1 + tl
        mem = max(v * 2.0 ** Bg, data)
        rows.append((max(data, time, mem), Bg, data, time, mem, B2, B3))
    rows.sort()
    obj, Bg, d, t, m, B2, B3 = rows[0]
    print("v=%d: Bg=%d data=2^%.3f time=2^%.3f mem=2^%.3f (B2=%d B3=%d)"
          % (v, Bg, math.log2(d), math.log2(t), math.log2(m), B2, B3))

# limit for v -> infinity, i.e. N1 -> 0
best = None
for Bg in range(1, L + 1):
    _, B2, B3, dl, tl = later_passes(L - Bg)
    time = Bg * 2.0 ** Bg + tl
    data = dl
    mem = max(2.0 ** Bg, data)  # memory at fixed v; it grows with v for the metric
    obj = max(data, time, mem)
    if best is None or obj < best[0]:
        best = (obj, Bg, data, time, mem)
obj, Bg, d, t, m = best
print("limit v->inf (time only): Bg=%d data=2^%.3f time=2^%.3f"
      % (Bg, math.log2(d), math.log2(t)))
print("interior minimum of the metric at v of order 2^%.2f (memory v 2^Bg reaches the time)"
      % (math.log2(t) - Bg))

# curve of the largest figure at v = 26 around its minimum (Appendix A)
S26 = 2.0 ** SEI[26]
print("v=26 curve:")
for Bg in range(210, 219):
    N1 = 4 * Bg * math.log(2) / S26
    T1 = 2 * N1 + Bg * 2.0 ** Bg
    _, B2, B3, dl, tl = later_passes(L - Bg)
    data = max(N1, dl)
    time = T1 + tl
    mem = max(26 * 2.0 ** Bg, data)
    print("%d & %.3f & %.3f & %.3f & %.3f" % (Bg, math.log2(data), math.log2(time), math.log2(mem), math.log2(max(data, time, mem))))
