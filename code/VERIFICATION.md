# Verification record

This file maps each computational statement of the paper to the script that
reproduces it. The record was produced from a clean copy of the repository on
2026-11-13 with Python 3.10.11 and numpy 2.0.2 on Windows, and with the .NET SDK
9.0.200 for the two C# enumerations. Every script below exited with status 0 and
reproduced the values shown. `code/README.md` gives the commands and the full
per-script table. The scripts added on 2026-11-18, `mask_search/correlation/triple_count.py`
and `attack_complexity/vectorial/tie_sensitivity.py`, and the scripts added on
2026-11-21 in `reduced_scale/`, were rerun on those dates with the same interpreter
(the reduced-scale scripts require only Python 3.10 and numpy).

## Mask search

| Paper object | Script | Observed output |
| --- | --- | --- |
| Table 1, distribution of the absolute values of the LAT and ILAT | `mask_search/correlation/verify.py` | 4590, 12240, 9180, 10200, 8670, 6120, 9180, 4080 and 1275 for the values 0, 4, 8, 12, 16, 20, 24, 28 and 32; 65535 entries, the (0,0) entry 256 excluded |
| Eq. (11), the optimal triple | `verify.py`, `check_triple.py`, `unit_tau.py` | tau = 2401/2^66 = 2^-54.770580312; the 18 absolute table values are 14 of 32 and 4 of 28 |
| The designers' trail | `verify.py` | tau = 1225/2^66 = 2^-55.741433966; the 18 absolute table values are 14 of 32, two of 28 and two of 20; mu(A)=mu(C)=mu(B L^-1)=0x0d0d and mu(B)=0x0505 |
| Ratio of the two correlations | `verify.py` | 49/25 = 2^0.970853654, that is 0.9709 bits |
| Section 3.3 support enumeration, Table 2 | `global_support_enumeration.cs`, `global_support_verify.cs` | eq18 all=158317600 screened=373544 feasible=4 same=4; lower all=314668232 screened=552044 feasible=0 |
| Exact maxima of the four supports | `global_optimum_scan.py` | 0x7070: -54.770580312, one pair with 100 triples; 0x0d0d, 0xd0d0 and 0x0707: -55.155870468 |
| Section 3.2 dimensions | `structure_check.py` | dim B = 16, dim C = 8, |Theta| = 2^24 for all four supports |
| Profile of the optimal support | `profile.py`, `verify_exact.py` | 11093111 nonzero pairs, 3295 distinct levels, range 2^-73.66 to 2^-54.771 |
| Head of the profile, Table 4 | `levels.py` | levels 0 to 7 (22 masks) plus the four masks at level 8 (-55.76355) |
| The 26-mask span V_26 | `sei_v26.py` | rank 26; SEI_2(V_26) = -216.797094; SEI_3(V_26) = -327.184547; log2 N at v=26 = 226.0100 |
| Section 4.1 pair counts and the approximation count of Section 4.2 | `triple_count.py` | nonzero pairs 11,093,111 / 11,193,038 / 11,032,069 / 11,082,491 (as in the paper); approximations 538,930,000 = 2^29.005523 (optimal support) and 2,154,559,625 = 2^31.004746 (four supports); the necessary condition of Theorem 4 holds |
| The earlier 24-mask span (auxiliary; the W of the paper is V_26) | `verify_list.py` | rank 24; 6138951 nonzero-correlation elements; log2 sum tau^4 = -216.8436670; log2 sum tau^6 = -327.2032269 |
| Section 4.3 masking kernels | `verify_mask_relations.py` | T=(0,2,3): rank 383, nullity 1; T=(0,1,2,4): rank 480, nullity 32 |

## Attack complexity

| Paper object | Script | Observed output |
| --- | --- | --- |
| Table 3, its first row and Eq. (12) | `attack_complexity/single/audit_attack_costs.py` | Bg=241: M = 2^228.0559, N = 2^250.0279, T = 2^250.9192; designers' trail under the rule of Section 3.4, Bg=243: N = 2^250.9745, T = 2^252.2014; designers' rule, Bg=248: N = 2^248.4866, T = 2^257.1863, M = 2^249.2637; optimum under the designers' rule, Bg=246: N = 2^247.5401, T = 2^255.8960, M = 2^247.9665 |
| Section 4.1: profiles of the four supports | `mask_search/correlation/profile_supports.py` | optimal 11093111 pairs / 3295 levels / -54.771; others 11193038/3165, 11032069/3266, 11082491/3299, all -55.156 |
| Section 4.2: capacities for v=24..29 (auxiliary) | `attack_complexity/vectorial/bg_scan_v.py`, `coset_v27.py`, `coset_v28.py` | auxiliary robustness scan with the corrected first-stage accounting, not quoted in the paper; SEI_2(V_29)=2^-216.726807 |
| Section 4.2: supports outside the profile | `attack_complexity/vectorial/coset_other.py`, `coset_all.py` | each further triple adds 0.098 bits to log2 SEI_2; all three add 0.2768 bits; the first-stage data falls by the same amount and `bg_opt.py` gives the change of the totals (time 0.167 bits lower, data and memory unchanged) |
| Section 4.2: tie-breaking rule of the 26 representatives | `attack_complexity/vectorial/tie_sensitivity.py` | 78 tied active bytes; rank 26 for every rule; log2 SEI_2 = -216.797094 (smallest maximizer, the stored basis), -216.798145 (largest maximizer), -216.829548 (random seed); log2 SEI_3 = -327.184547 / -327.184618 / -327.193051; the change of log2 N is at most 0.017 bits |
| Single-attack parity-count sensitivity (footnote 3 in Section 3.4) | hand calculation from the formulas of Section 3.4 with the sample-generation coefficient tripled | fixed \(B_g\): total time \(2^{252.255}\) and \(2^{253.344}\), improvement 1.089 bits; reoptimized \(B_g=242\) and \(244\): \(2^{251.986}\) and \(2^{253.255}\), improvement 1.269 bits |
| Operation-count sensitivity (footnote 4 in Section 4.2) | hand calculation from the formulas of Section 4.2 with three parities per clock | 156N sample generation: \(B_g=75\), \(T_1=2^{231.801}\), total time \(2^{232.109}\), data \(2^{228.938}\), gains 18.810 bits in time and 21.090 in data |
| Table 1 row (d), Table 5, Section 4.2 final vectorial costs | `attack_complexity/vectorial/cost_v29.py` | v = 29, B_g = 239; data `2^247.455`, time `2^248.522`, memory `2^247.455`; gains 2.404/2.288/2.404; overall success probability 0.89998 |
| Table 5, Table 6, Sections 4.2 and 4.3 | `attack_complexity/vectorial/cost.py`, `bg_opt.py` and `sei_v26.py` | v=1: 2^-219.082; v=26: SEI_2 = -216.797094, SEI_3 = -327.184547, log2 N at Bg=214 = 226.010; bg_opt.py: B_g = 78, first stage N = 224.554, T = 230.309, M = 82.700; later passes 216 bits (N = 228.927, T = 229.538) and 218 bits (N = 118.396, T = 224.907), total data 228.927, time 230.996, memory 228.927, gain 19.923 bits in time and 21.101 in data over the single approximation; with SEI_2 = 2^-216.520327 the time falls 0.167 bits and data and memory are unchanged |

## Supporting checks

`ties.py`, `fourth2.py`, `t24.py`, `verify_W2.py`, `screen24.py`, `verify_space.py`,
`search_space.py`, `pair_search.py`,
`push_k24.py`, `rank_check.py`, `crosscheck_sei.py`, `feas.py` and `feas2.py`
recompute internal quantities that the argument consumes (bytewise maximisers,
power sums, basis selections and feasibility witnesses). Their expected values
are listed in `code/README.md`.

## Reduced-scale validation (supplementary)

| Paper object | Script | Observed output |
| --- | --- | --- |
| Supplementary reduced-scale checks, structure of the analogue | `reduced_scale/verify.py` | all checks V1 to V9 pass; S-box and multipliers as specified; the two forms of e(t) agree bitwise; the piling-up identity holds exactly at clock zero by exhaustive enumeration over 2^24 register states; the paper-convention formula checks by exhaustive enumeration |
| Supplementary reduced-scale checks, mask search | `reduced_scale/search_masks.py` | 111 triples with log2 abs(tau) = -9.000; three linearly independent basis triples; predicted span capacity 2^-34.415 |
| Table 8, single-clock estimates | `reduced_scale/monte_carlo.py` | at N=2^26: t1 +1.9688e-3, t2 +1.9155e-3, t3 +1.8990e-3, w5 +1.016e-4, w7 -5.52e-5, standard error 2^-13; every estimate within one standard error of the prediction |
| Table 8, folded correlations and capacity | `reduced_scale/exact_folded.py` | exact over 2^24 register states: d=1 ratios -4.0, -6.0, -50.0, -409.6 and 0; d=3 ratio exactly one; capacity d=1 2^-24.681 against the prediction 2^-34.415, and d=3 2^-34.415 |

## Reproducibility notes

- `data/w24.log` is written without the run time, so it is byte-identical
  across runs of `verify_list.py`. The printed summary still reports the time.
- `code/mask_search/search/sei_nested.json` and `sei_nested_24.json` are
  regenerated byte-identically by `verify_space.py` and `push_k24.py`.
- The two C# files define the static methods `SnowSupports.Run()` and
  `IndependentEnumeration.Run()` and no entry point, so each is compiled with a
  one-line `Main`, as described in `code/README.md`.
- An independent implementation written for this verification, and not part of
  this repository, reproduced the ten mask-search checks above on 2026-11-13,
  including the exact rationals 2401/2^66, 1225/2^66 and 49/25.
