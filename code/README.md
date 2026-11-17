# Code for "Linear Cryptanalysis of SNOW-Axn"

The code reproduces the computational results of the paper "Linear Cryptanalysis of SNOW-Axn". It is organised in two categories: the mask search,
which establishes the correlation results on the three-word approximation, and the attack
complexity, which gives the cost of the single-approximation, vectorial and distinguishing attacks.
Every script listed below was rerun from a clean copy on 2026-11-13, and the observed values are
recorded in `VERIFICATION.md`.

    code/
      README.md
      mask_search/
        correlation/      LAT and ILAT tables, the per-byte factorisation, and the
                          exact correlation of the optimal and designers' triples
        search/           support enumeration, exact optimum, feasibility and dimension
                          checks, and the search for the 24-mask space W
      attack_complexity/
        single/           single-approximation fast correlation attack
        vectorial/        vectorial fast correlation attack
        distinguishing/   distinguishing attack
      dev/                development history, not needed for reproduction

## Requirements and how to run

* Python 3.10 or newer with numpy. The last full reproduction used Python 3.10.11 and numpy 2.0.2 on Windows.
* The .NET SDK (or another C# compiler) for the two enumerations. Each file defines a static `Run()`
  method and no entry point, so compile it with a one-line `Main`. With the .NET SDK, for example,

        dotnet new console -o enum
        cp <path>/global_support_enumeration.cs enum/
        printf 'public static class Program { public static void Main() { SnowSupports.Run(); } }' > enum/Program.cs
        dotnet run --project enum

  For `global_support_verify.cs` the call is `IndependentEnumeration.Run()`.
* Any script runs from any directory, for example `python code/mask_search/correlation/verify.py`.
  Every script outside `dev/` locates `verify.py` and the data from its own path.

## Mask search

### correlation

`verify.py` rebuilds GF(2^8), the AES S-box, the linear layers L and sigma, and the LAT
and ILAT tables, factorises the three-block correlation per byte, and evaluates the
optimal triple of Eq. (11) and the designers' triple of Appendix A of the design
document. The designers' trail is reproduced box by box: 18 active S-boxes with table
absolute values 32, 28 and 20 and total correlation 2^-55.7414.

| Script | Paper object | Expected output |
| --- | --- | --- |
| `verify.py` | Table 1 (LAT and ILAT absolute-value distribution) and Eq. (11) (both triples) | `log2 tau* = -54.7706`, designers `-55.7414`, gap `0.9709` bits |
| `check_triple.py` | Eq. (11), recomputed from the mask bytes printed in the paper | `log2 tau=-54.770580` |
| `unit_tau.py` | Eq. (11), byte-by-byte unit check | `log2 tau = -54.77058031176958` |
| `verify_exact.py` | Section 4.1 correlation profile | 11,093,111 non-zero, max `log2\|tau\| = -54.7706` at one pair `(B,C)`, 3,295 levels |
| `profile.py` | Section 4.1 by exhaustive enumeration of Theta_mu | 11,093,111 non-zero, max `log2\|tau\| = -54.7706`, 3,295 distinct levels |
| `levels.py` | Table 4, the level list near the head | `all levels 0..7 (22 masks) + 2 of the 4 masks at level 8`; `uniquely determined: False` |
| `ties.py` | number of bytewise maximizers among the 24 strongest | 11,476 of 65,536 pairs have more than one maximiser; 72 tied active bytes among the 24 masks |
| `fourth2.py` | profile sums and the two cost rules behind Section 4.2 | `sum tau^2 = 2^-97.938`, `sum tau^4 = 2^-213.382`, `sum tau^6 = 2^-326.187` |
| `t24.py` | tau^2, tau^4 and the cross-term bound over the 24 strongest | `log2 sum tau^2 = -106.272`, `log2 sum tau^4 = -216.881`, relative bound `-69.003` |
| `verify_W2.py` | span of one 22 + 2 selection of the 24 strongest (used by `ties.py`) | `python verify_W2.py 22,24`: rank 24, 6,254,429 non-zero; the canonical (22,23) selection gives 6,138,951 and is enumerated by `verify_list.py` |

### search

`global_support_enumeration.cs` and `global_support_verify.cs` enumerate every support with
exactly eighteen active S-boxes and screen those with fewer, and both find the same four,
each feasible with activity 18.
`global_optimum_scan.py` evaluates the exact optimum of those four supports.
`feas.py`, `feas2.py` and `structure_check.py` check feasibility and the dimensions 16
and 8, and `screen24.py` screens the choice of the last two masks. The remaining scripts
search the nested spans of the profile and return the 24 triples that span W;
`verify_list.py` in `attack_complexity/distinguishing` enumerates the resulting span.

| Script | Paper object | Expected output |
| --- | --- | --- |
| `global_support_enumeration.cs` | Section 3.3 support enumeration, Table 2 | `eq18 all=158317600 screened=373544 feasible=4 same=4`, `lower all=314668232 screened=552044 feasible=0` |
| `global_support_verify.cs` | the same enumeration by row-space membership | the same totals and the same four hits |
| `global_support_matrix.py` | byte matrices and the four MDS blocks behind the proof | `256` basis pullbacks agree, all `276` minors nonzero |
| `global_optimum_scan.py` | Section 3.3 and Table 2, exact optimum of the four supports | `mu=7070 ... log2tau=-54.770580312`, one optimal pair with 100 triples; three at `-55.155870468` |
| `feas.py` | random sparse feasibility of the 0/1 model | `feasible witness found: False` |
| `feas2.py` | structured feasibility with gamma = sigma(L(alpha)) | `structured hits: 6` |
| `structure_check.py` | Section 3.2 dimensions of B_mu and C_mu | `dim B = 16, dim C = 8, \|Theta\| = 2^24` for all four patterns |
| `screen24.py` | rank screening of the six 22 + 2 selections | `rank384 = 24` for all six, `rank(B,C) = 20` or `21` |
| `verify_space.py` | the two basis-selection strategies; writes `sei_nested.json` next to itself | `nested best: k=1 log2 T=221.5536`, `greedy best: k=1 log2 T=221.5536` |
| `search_space.py` | nested spans of the k strongest masks | `best k=1 log2 T=221.5536 log2 SEI=-219.0823` |
| `pair_search.py` | pairs and triples from a 400-mask pool | `best pair: SEI -218.6562`; `best triple: SEI -218.3277` |
| `push_k24.py` | extends the nested search to k = 24; writes `sei_nested_24.json` next to itself | `k=21 ... k=24 log2 SEI = -216.9048 ... -216.8437` |
| `rank_check.py` | GF(2) rank of the 24 strongest complete triples | `rank of the 24 strongest masks = 24` |
| `crosscheck_sei.py` | independent SEI for k = 1, 2, 3 and sampled k = 24 | `k=1 direct SEI = -219.082321`, `k=2 = -218.656208`, `k=3 = -218.327675` |

## Attack complexity

Each subsection names the attack and its cost, its references and the paper section. All
costs are base-2 logarithms of data N, time T and memory.

### single

Single-approximation fast correlation attack, Meier and Staffelbach (1989) and Selcuk (2008).
Gives the cost of Section 3.4.

| Script | Computes | Expected output |
| --- | --- | --- |
| `audit_attack_costs.py` | Table 3, its first row under the designers' rule, and the statistical constants of the first stage | `bg=241 ... time=250.919215250`; designers' trail under the rule of Section 3.4, `bg=243 ... time=252.201402061`; designers' rule, `bg=248 ... time=257.186283588`; optimum under the designers' rule, `bg=246 ... time=255.896018098` |

### vectorial

Vectorial fast correlation attack, Zhang, Liu, Meier, Sun, Feng and Wu (ASIACRYPT
2025). Gives the cost of Section 4.2. The statistic is a vector over the space, one sample carries
the whole vector, and the rule assigns N = 4 kappa ln2 / SEI(V_v) keystream words,
2N + kappa 2^kappa in time and v 2^kappa words of memory for the transform, where kappa
is the number of state bits decoded in one step and v is the dimension of the space. The working point kappa = 214 is optimized
in bg_opt.py, and the total memory is the stored samples.

| Script | Computes | Expected output |
| --- | --- | --- |
| `cost.py` | Table 5 (scan for v = 1..26 at the working point) | `v = 26: log2 N = 226.010` |
| `bg_opt.py` | the B_g scan of Section 4.2 with the two later passes | `B_g = 214`, total data `227.857`, time `228.918`, memory `227.857`, gain `22.001` bits in time, `22.171` in data |
| `coset_other.py` | per-support effect on the capacity of W (Section 4.2) | each strongest further triple gives `coset sum4 = 2^-220.623482`; `log2 SEI_2` rises by `0.098` |
| `bg_scan_v.py` | auxiliary robustness scan over v and B_g with the two later passes, not quoted in the paper | `v = 24..29: Bg = 214, data = 227.857`; time `228.931/228.924/228.918/228.911/228.903/228.900`; time limit `2^228.471`; curve of the largest figure at v=26 around `Bg = 214` (Appendix A, Table 7) |
| `coset_v28.py` | capacity of V_28 (Section 4.2) | `coset sum4 = 2^-222.395499`; `SEI_2(V_28) = 2^-216.741826`; first-stage data falls 0.029 bits |
| `coset_v27.py` | capacity of V_27 (Section 4.2) | `coset sum4 = 2^-222.560324`; `SEI_2(V_27) = 2^-216.770773`; first-stage data falls 0.026 bits |
| `coset_all.py` | total effect of the three further feasible supports on the capacity of W (Section 4.2) | singles `2^-220.623482`, pairs `2^-441.246964`/`2^-373.606336`, triple `2^-594.229818`; total `log2 SEI_2` rise `0.277`; data fall `0.277` bits |
| `sei_v26.py` | Table 5 rows v = 25, 26 and the numbers of Sections 4.2 and 4.3 (power sums of V_26) | `SEI_2(V_26) = 2^-216.797094`, `SEI_3(V_26) = 2^-327.184547`, `log2 N = 226.010` |
| `sei_nested.json`, `sei_nested_24.json` | cached capacities for v = 1..20 and v = 21..24, read by `cost.py` | the same files that `verify_space.py` and `push_k24.py` regenerate |

### distinguishing

Distinguishing attack, Nyberg and Hermelin (2007) and Yang, Johansson and Maximov (2020).
Gives the cost of Section 4.3.

| Script | Computes | Expected output |
| --- | --- | --- |
| `profile_supports.py` | profiles of the four feasible supports (Section 4.1) | optimal `11093111/3295/-54.771`; others `11193038/3165`, `11032069/3266`, `11082491/3299`, all `-55.156` |
| `verify_list.py` | enumerates the earlier 24-mask span of `../data/w24.json` (auxiliary; the W of the paper is V_26 in `sei_v26.py`) | rank 24, 6,138,951 non-zero, `log2_sum4=-216.8436669675`, `log2_sum6=-327.2032268545` |
| `verify_mask_relations.py` | separated-time linear masking kernel ranks | `T=(0, 2, 3) rank=383 nullity=1`, `T=(0, 1, 2, 4) rank=480 nullity=32` |

## Data

`../data/w24.json` holds the 24 triples of the earlier 24-mask space (auxiliary; the W of the paper is V_26), with byte 0 at the left end and the
selection rule in the file. Its SHA-256 is
`68ef228943f9f4a2fafa571e28ce43bc06a9437b68662c12b22b8fe927a256e1`, and `../data/w24.log`
is the output of the last run of `verify_list.py`. The run time is printed but not written,
so the log is byte-identical across runs.

`../data/w26.json` holds the 26 triples that span V_26, the 22 of levels 0 to 7 and the
four of level 8, with the selection rule in the file. Its SHA-256 is
`112361d5c62f0104bda51a64240812d47f597032b73e13fdcd6ba5a497f4bbeb`; `sei_v26.py`
enumerates the 26-dimensional span in four cosets of dimension 24.

The cached capacities are `attack_complexity/vectorial/sei_nested.json` for k = 1..20 and
`sei_nested_24.json` for k = 21..24. `verify_space.py` and `push_k24.py` regenerate the
same files, byte for byte, next to themselves in `mask_search/search/`.

## Runtime

All Python scripts require Python 3.10 or newer with numpy (the last full run used Python
3.10.11 and numpy 2.0.2). Every script outside `dev/` locates `verify.py` and the data
from its own path, so it can be run from any directory. `verify_list.py` enumerates 2^24
sums and takes about two minutes, `push_k24.py` about one minute, and `verify_exact.py`,
`profile.py`, `t24.py`, `screen24.py` and the nested searches tens of seconds. The two C#
enumerations take about ten seconds each including the build, and the other scripts
finish in seconds.

## Licence

The code is released under the MIT licence; see `../LICENSE`. The paper text and figures
are released under CC-BY-4.0.

`dev/` holds development history, one-off probes and superseded drafts; it is not needed
for reproduction. It also keeps `vectorial_snow.py` and `verify_synchronous_collision.py`,
which estimate and check the earlier paired-clock model that the paper no longer uses.
