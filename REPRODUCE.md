# Minimal reproduction of the current paper numbers

Environment: Python 3.10+ with numpy.

## 1. Exact three-clock optimum

    python code/mask_search/correlation/verify.py
    python code/mask_search/search/global_optimum_scan.py

Expected: `log2 tau* = -54.770580312`, designers `-55.741...`, four feasible
supports, and the four support optima `-54.770580312` / `-55.155870468`.

## 2. Correlation attack cost

    python code/attack_complexity/single/audit_attack_costs.py

Expected: the single-attack figures used in Table 1.

## 3. Final vectorial attack (v = 29)

    python code/attack_complexity/vectorial/cost_v29.py

Expected:

    basis triples: 29
    SEI_2(W) = 2^-211.540788
    optimal B_g = 239
    data  = 2^247.456
    time  = 2^248.523
    memory= 2^247.456
    gains: data 2.403, time 2.287, memory 2.403

`data/w29_opt.json` is the basis. The v = 29 point is the best point of the
joint (v, B_g) scan, not a global optimality claim.

## 4. Distinguisher

    python code/attack_complexity/vectorial/sei_v26.py

Expected: `SEI_2(V_26) = 2^-216.797094`, `SEI_3(V_26) = 2^-327.184547`.
