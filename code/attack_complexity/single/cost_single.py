#!/usr/bin/env python3
"""Cost of the single-approximation attack (Table 1 row (c))."""
import json
import math

L = 512
PS = 29.0 / 30.0
TAU = 2401.0 / 2.0 ** 66
BG = 241
T3 = 2.0 ** 228.560


def qinv(p):
    lo, hi = 0.0, 60.0
    for _ in range(200):
        mid = (lo + hi) / 2.0
        if 0.5 * math.erfc(mid / math.sqrt(2.0)) > p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def main():
    rho = TAU * TAU
    M = (qinv(1 - PS) + qinv(2.0 ** (-BG - 1))) ** 2 / rho ** 2
    N = math.sqrt(M * 2.0 ** (L + 1 - BG))
    T1 = 2.0 * math.sqrt(M * 2.0 ** (L - BG))
    T2 = BG * 2.0 ** BG + M
    total = T1 + T2 + T3
    print("M     = 2^%.6f" % math.log2(M))
    print("N     = 2^%.6f" % math.log2(N))
    print("T1    = 2^%.6f" % math.log2(T1))
    print("T2    = 2^%.6f" % math.log2(T2))
    print("T3    = 2^%.6f" % math.log2(T3))
    print("T_tot = 2^%.6f" % math.log2(total))
    assert abs(math.log2(N) - 249.859) < 1e-3
    assert abs(math.log2(total) - 250.810) < 1e-3


if __name__ == "__main__":
    main()
