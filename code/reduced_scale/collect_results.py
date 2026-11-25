#!/usr/bin/env python3.10
"""Assemble the report tables from the JSON outputs of the experiment.

Reads masks.json, exact_results.json, mc_n24.json, mc_n26.json (all produced
by the scripts in this directory) and prints Markdown tables.

Run:  python3.10 collect_results.py
"""

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    with open(os.path.join(HERE, name), encoding="utf-8") as f:
        return json.load(f)


def log2(x):
    return float(np.log2(abs(x))) if x else float("-inf")


def main():
    masks = load("masks.json")
    exact = load("exact_results.json")
    m24 = load("mc_n24.json")
    m26 = load("mc_n26.json")
    basis = masks["basis"]
    combos = masks["span"]

    print("## T1 masks (top by |tau|)\n")
    print("| label | a | b | c | c'=sigma^-1(c) | tau_num | tau | log2\\|tau\\| |")
    print("|---|---|---|---|---|---|---|---|")
    labels = {1: "t1", 2: "t2", 3: "t1^t2", 4: "t3", 5: "t1^t3", 6: "t2^t3",
              7: "t1^t2^t3"}
    for i, r in enumerate(masks["top"][:8]):
        tag = ""
        for bi, br in enumerate(basis):
            if (r["a"], r["b"], r["c"]) == (br["a"], br["b"], br["c"]):
                tag = " **%s**" % labels[1 << bi] if bi < 3 else ""
        print("| top%d%s | 0x%02x | 0x%02x | 0x%02x | 0x%02x | %d | %+0.9f | %.6f |"
              % (i + 1, tag, r["a"], r["b"], r["c"], r["cprime"], r["tau_num"],
                 r["tau"], log2(r["tau"])))
    print()
    print("Basis (linearly independent triples):")
    for bi in range(3):
        r = basis[bi]
        print("- t%d = (0x%02x, 0x%02x, 0x%02x), tau = %+0.9f = 2^%.3f"
              % (bi + 1, r["a"], r["b"], r["c"], r["tau"], log2(r["tau"])))

    print()
    print("## T2 single-clock correlations\n")
    print("| label | a,b,c | tau_pred (log2) | t=0 exact | t=0 MC | t=1..6 MC |"
          " t=1..6 SE | dev/SE |")
    print("|---|---|---|---|---|---|---|---|")
    for i, (r24, rex) in enumerate(zip(m24["combos"], exact["combos"])):
        tau = r24["tau_pred"]
        se = r24["single_t1_6_se"]
        dev = (r24["single_t1_6"] - tau) / se if se else 0.0
        print("| %s | 0x%02x,0x%02x,0x%02x | %.6f | %+0.9f | %+0.9f | %+0.9f |"
              " %.3e | %+.2f |"
              % (labels[r24["w"]], r24["a"], r24["b"], r24["c"], log2(tau),
                 rex["tau_exact"], r24["single_t0"], r24["single_t1_6"], se, dev))

    print()
    print("## T3 folded correlations\n")
    print("| label | d | tau^2 pred (log2) | C2 exact | log2\\|C2\\| |"
          " MC 2^24 | SE | MC 2^26 | SE | exact/pred | bits |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    keys = {1: "C2_1", 2: "C2_2", 3: "C2_3"}
    for i, (rex, r24, r26) in enumerate(zip(exact["combos"], m24["combos"],
                                            m26["combos"])):
        for d in (1, 2, 3):
            ex = rex[keys[d]]
            f24 = r24["folded"][d - 1]
            f26 = r26["folded"][d - 1]
            tau2 = rex["tau_pred"] ** 2
            ratio = ex / tau2 if tau2 else float("nan")
            bits = log2(ex) - log2(tau2) if (ex and tau2) else float("nan")
            print("| %s | %d | %.6f | %+0.6e | %.3f | %+0.6e | %.2e | %+0.6e |"
                  " %.2e | %s | %s |"
                  % (labels[rex["w"]], d, log2(tau2), ex, log2(ex),
                     f24["C2"], f24["se"], f26["C2"], f26["se"],
                     ("%.3f" % ratio) if tau2 else "n/a",
                     ("%+.3f" % bits) if (ex and tau2) else "n/a"))

    print()
    print("## T4 capacity SEI_2 (v=3 span)\n")
    print("| quantity | value | log2 |")
    print("|---|---|---|")
    print("| predicted sum_w tau^4 | %.6e | %.6f |"
          % (exact["capacity"]["SEI2_pred"], log2(exact["capacity"]["SEI2_pred"])))
    print("| exact, folded offset d=1 | %.6e | %.6f |"
          % (exact["capacity"]["SEI2_exact"], log2(exact["capacity"]["SEI2_exact"])))
    print("| exact, folded offset d=3 | %.6e | %.6f |"
          % (exact["capacity"]["SEI2_d3"], log2(exact["capacity"]["SEI2_d3"])))
    print("| MC 2^24 (d=1, noise-dominated) | %.6e | %.6f |"
          % (m24["capacity"]["empirical_SEI2"], log2(m24["capacity"]["empirical_SEI2"])))
    print("| MC 2^26 (d=1, noise-dominated) | %.6e | %.6f |"
          % (m26["capacity"]["empirical_SEI2"], log2(m26["capacity"]["empirical_SEI2"])))

    print()
    print("## T5 different-clock independence test (MC 2^26)\n")
    print("| label | d | C2 empirical | SE | product of single-clock means |"
          " difference | SE | diff/SE |")
    print("|---|---|---|---|---|---|---|---|")
    for r26 in m26["combos"]:
        for ind in r26["independence"]:
            print("| %s | %d | %+0.9f | %.2e | %+0.3e | %+0.9f | %.2e | %+.2f |"
                  % (labels[r26["w"]], ind["d"], ind["C2"], ind["C2_se"],
                     ind["single_product"], ind["difference"],
                     ind["difference_se"], ind["diff_over_se"]))

    print()
    print("## T6 exact deviation of the folded check from tau^2\n")
    print("| label | tau | d=1: C2/tau^2 | d=1 bits | d=2: C2/tau^2 |"
          " d=2 bits | d>=3 |")
    print("|---|---|---|---|---|---|---|")
    for rex in exact["combos"]:
        tau = rex["tau_pred"]
        e1, e2 = rex["C2_1"], rex["C2_2"]
        r1 = e1 / tau ** 2 if tau else float("nan")
        r2 = e2 / tau ** 2 if tau else float("nan")
        b1 = log2(e1) - 2 * log2(tau) if (e1 and tau) else float("nan")
        b2 = log2(e2) - 2 * log2(tau) if (e2 and tau) else float("nan")
        print("| %s | %+0.6e | %s | %s | %s | %s | C2=tau^2 |"
              % (labels[rex["w"]], tau,
                 ("%.3f" % r1) if tau else "n/a",
                 ("%+.3f" % b1) if (e1 and tau) else ("n/a" if not e1 else ""),
                 ("%.3f" % r2) if tau else "n/a",
                 ("%+.3f" % b2) if (e2 and tau) else "n/a"))


if __name__ == "__main__":
    main()
