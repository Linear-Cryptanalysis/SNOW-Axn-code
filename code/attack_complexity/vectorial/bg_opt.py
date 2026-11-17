"""Optimized working point of the vectorized attack at v = 26.

Computes the cost rule of Theorem 4 of Zhang et al. (ASIACRYPT 2025) for every
B_g, with the two later passes of Section 3.4 included, and prints the optimum
and the value at B_g = 240.

Rule: N = 4 B_g ln 2 / SEI_2(V_v), T = 2N + B_g 2^B_g, M = v 2^B_g.

Inputs: none. The capacity of V_26 is 2^-216.797094.
Outputs: the optimum B_g with the data, time and memory of the full attack, and
         the value at B_g = 240.
Paper object: Section 4.2 and the figures quoted in the abstract.
"""

import math

PS = 0.999992
L = 512
TAU = 2401 / 2 ** 66
SEI2 = 2.0 ** -216.797094
V = 26                       # dimension of the space W = V_26
SINGLE_TIME = 250.919
SINGLE_DATA = 250.028


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


def Meq(c, B):
    return (1 / c ** 2) * (q1 + qinv(2.0 ** (-B - 1))) ** 2


def dec(B, M):
    if B <= math.log2(M):
        return B * 2.0 ** B + M
    return 2.0 ** B * (1 + math.ceil(math.log2(M)))


def later_passes(rem):
    best = None
    for B2 in range(rem + 1):
        B3 = rem - B2
        M2 = Meq(TAU * TAU, B2)
        N2 = 2 * M2
        T2 = N2 + dec(B2, M2)
        if B3 > 0:
            M3 = Meq(TAU, B3)
            N3 = M3
            T3 = N3 + dec(B3, M3)
        else:
            N3, T3 = 0.0, 0.0
        data = max(N2, N3)
        time = T2 + T3
        if best is None or max(time, data) < best[0]:
            best = (max(time, data), B2, B3, data, time)
    return best


def cost(Bg):
    N1 = 4 * Bg * math.log(2) / SEI2
    T1 = 2 * N1 + Bg * 2.0 ** Bg
    M1 = V * 2.0 ** Bg
    _, B2, B3, dl, tl = later_passes(L - Bg)
    data = max(N1, dl)
    time = T1 + tl
    mem = max(M1, data)
    return (max(data, time, mem), data, time, mem, N1, T1, M1, B2, B3)


rows = [(cost(Bg)[0], Bg) for Bg in range(1, L + 1)]
rows.sort()
obj, Bg = rows[0]
obj, data, time, mem, N1, T1, M1, B2, B3 = cost(Bg)
print("optimum B_g = %d" % Bg)
print("  first stage: log2 N = %.3f  log2 T = %.3f  log2 M = %.3f" % (math.log2(N1), math.log2(T1), math.log2(M1)))
print("  two later passes: B2 = %d  B3 = %d" % (B2, B3))
print("  totals: log2 data = %.3f  log2 time = %.3f  log2 memory = %.3f" % (math.log2(data), math.log2(time), math.log2(mem)))
print("  gain over the single attack: %.3f bits in time, %.3f bits in data"
      % (SINGLE_TIME - math.log2(time), SINGLE_DATA - math.log2(data)))
obj240, data240, time240, mem240 = cost(240)[:4]
print("B_g = 240: log2 data = %.3f  log2 time = %.3f  log2 memory = %.3f"
      % (math.log2(data240), math.log2(time240), math.log2(mem240)))
