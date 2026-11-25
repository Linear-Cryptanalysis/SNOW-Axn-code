"""Estimate the cost of the vectorial attack of Section 4.2.

Convention of the paper. The statistic is a vector over the space V_v, so one
sample carries the whole vector and no pair of clocks has to collide. The rule
of Zhang, Liu, Meier, Sun, Feng and Wu (ASIACRYPT 2025) then assigns

    N = 4 kappa ln 2 / SEI(V_v),   T = 54N + kappa 2^kappa,   memory = v 2^kappa,

where kappa is the number of state bits decoded in one step and v is the
dimension of the space, which enters through the capacity and through the memory. The working
point kappa = 214 and the two later passes are optimized in bg_opt.py.

Inputs: sei_nested.json (v = 1..20), sei_nested_24.json (v = 21..24) and the
        values for v = 25, 26 printed by sei_v26.py.
Outputs: the dimension scan of Table 5 at the working point kappa = 214.
Paper object: Table 5 and the equations of Section 4.2.
"""

import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
KAPPA = 214                      # state bits decoded in one step, optimized in bg_opt.py
TAIL = [-216.820106, -216.797094]  # v = 25, 26


def log2_data(log2sei):
    """4 kappa ln 2 / SEI, in base-two logarithm."""
    return math.log2(4.0 * KAPPA * math.log(2.0)) - log2sei


def main():
    sei = (json.load(open(os.path.join(HERE, 'sei_nested.json'))) +
           json.load(open(os.path.join(HERE, 'sei_nested_24.json'))) + TAIL)
    assert len(sei) == 26
    print('%2s %12s %10s' % ('v', 'log2 SEI', 'log2 N'))
    best_v, best_n = 0, None
    for v, s in enumerate(sei, start=1):
        n = log2_data(s)
        print('%2d %12.3f %10.3f' % (v, s, n))
        if best_n is None or n < best_n:
            best_v, best_n = v, n
    print()
    print('kappa = %d  ->  log2 mem of the transform = %.3f (v = 26)' % (KAPPA, math.log2(26) + KAPPA))
    print('v = 26: log2 N = %.3f' % log2_data(sei[-1]))
    print()
    print('paper, Table 5 column at kappa = 214: v = 26 gives 226.010')
    assert best_v == 26
    assert abs(log2_data(sei[-1]) - 226.010) < 5e-4
    print('checks passed')


if __name__ == '__main__':
    main()
