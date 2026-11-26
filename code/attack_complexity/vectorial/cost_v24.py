#!/usr/bin/env python3
"""Cost of the 24-dimensional optimal vectorial attack (Section 4.2).

Input: data/w24_opt.json (deterministic optimal basis, rank 24).
Output: SEI_2(W), optimal B_g, data/time/memory, and gains over the single attack.

The two later passes use p_s = 0.999992; their cost is subdominant.
"""
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
DATA = os.path.join(ROOT, "data", "w24_opt.json")

L = 512
PS1 = 0.9
PS2 = 0.999992
SINGLE = (249.859, 250.810, 249.859)
LATER_DATA = 2.0 ** 226.358
LATER_TIME = 2.0 ** 229.331


def qinv(p):
    lo, hi = 0.0, 60.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if 0.5 * math.erfc(mid / math.sqrt(2.0)) > p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def main():
    meta = json.load(open(DATA, encoding="utf-8"))
    sei = 2.0 ** meta["log2_SEI2"]
    q1 = qinv(1 - PS1)
    best = None
    for bg in range(1, L + 1):
        M = 4 * bg * math.log(2) / sei
        nc = math.sqrt(M * 2.0 ** (L + 1 - bg))
        t1 = 2 * math.sqrt(M * 2.0 ** (L - bg)) + bg * 2.0 ** bg
        data = max(nc, LATER_DATA)
        time = t1 + LATER_TIME
        mem = max(nc, 24 * 2.0 ** bg)
        obj = max(data, time, mem)
        if best is None or obj < best[0]:
            best = (obj, bg, data, time, mem, M, nc, t1)
    _, bg, data, time, mem, M, nc, t1 = best
    print("basis triples: 24")
    print("SEI_2(W) = 2^%.6f" % math.log2(sei))
    print("optimal B_g = %d" % bg)
    print("M     = 2^%.3f" % math.log2(M))
    print("N_c   = 2^%.3f" % math.log2(nc))
    print("T_1   = 2^%.3f" % math.log2(t1))
    print("data  = 2^%.3f" % math.log2(data))
    print("time  = 2^%.3f" % math.log2(time))
    print("memory= 2^%.3f" % math.log2(mem))
    print("gains: data %.3f, time %.3f, memory %.3f" % (
        SINGLE[0] - math.log2(data),
        SINGLE[1] - math.log2(time),
        SINGLE[2] - math.log2(mem),
    ))


if __name__ == "__main__":
    main()
