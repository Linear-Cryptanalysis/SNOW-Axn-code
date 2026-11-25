# -*- coding: utf-8 -*-
"""Cost of the vectorial attack on the 29-dimensional space W (paper Section 4.2).

Model of the paper:
  M    = 4 B_g ln 2 / SEI_2(W)
  N_c  = sqrt(M 2^(l+1-B_g))
  T_1  = 2 (M 2^(l-B_g))^(1/2) + B_g 2^B_g
  M_1  = max(N_c, v 2^B_g)
The two later passes use p_s = 0.999992 and contribute data 2^227.356
and time 2^228.515; the total data/time are dominated by the first stage.
"""
import math

L = 512
V = 29
SEI2 = 2.0 ** -216.520327
LATER_DATA = 2.0 ** 227.356
LATER_TIME = 2.0 ** 228.515
SINGLE_DATA = 249.859
SINGLE_TIME = 250.810

best = None
for Bg in range(1, L + 1):
    M = 4 * Bg * math.log(2) / SEI2
    Nc = math.sqrt(M * 2.0 ** (L + 1 - Bg))
    T1 = 2 * math.sqrt(M * 2.0 ** (L - Bg)) + Bg * 2.0 ** Bg
    mem = max(Nc, V * 2.0 ** Bg)
    data = max(Nc, LATER_DATA)
    time = T1 + LATER_TIME
    obj = max(data, time, mem)
    if best is None or obj < best[0]:
        best = (obj, Bg, data, time, mem, M, Nc, T1)

_, Bg, data, time, mem, M, Nc, T1 = best
print("W basis: v = %d, SEI_2 = 2^%.6f" % (V, math.log2(SEI2)))
print("optimum B_g = %d" % Bg)
for name, x in [("M", M), ("N_c", Nc), ("T_1", T1), ("data", data), ("time", time), ("memory", mem)]:
    print("%-8s 2^%.3f" % (name, math.log2(x)))
print("gains over single: data %.3f, time %.3f, memory %.3f"
      % (SINGLE_DATA - math.log2(data), SINGLE_TIME - math.log2(time),
         SINGLE_DATA - math.log2(mem)))
