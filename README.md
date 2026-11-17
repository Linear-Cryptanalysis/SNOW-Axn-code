# Code and data for "Linear Cryptanalysis of SNOW-Axn"

This repository holds the code and the data that reproduce the computational results of
the paper. The paper itself is not part of the repository.

SNOW-Axn is a 256-bit AEAD stream cipher of the SNOW family. Its designers assessed
resistance to linear cryptanalysis with a three-word approximation of correlation
2^-55.741, obtained from a heuristic search, which left open whether a better
approximation exists. The paper settles that model exactly: enumerating every support
with eighteen active S-boxes leaves four supports, and screening those with fewer
leaves none. Their exact evaluation gives the optimum 2^-54.771, which is 0.971 bits
above the designers' trail and 0.771 bits short of the 2^-54 bound implied by their count
of active S-boxes. With a known LFSR state the optimal approximation recovers the 512-bit
state at 2^250.919 in time and 2^250.028 in data and memory, while without a known state
the strongest approximations give a distinguisher on the 26-dimensional space W that needs
2^327.184 keystream words.

## Repository layout

    code/                              Python 3 and C# scripts; code/README.md maps each script to the
                                       paper object it reproduces
    code/mask_search/correlation/      LAT and ILAT tables, the per-byte factorisation,
                                       and the exact correlations of the two triples
    code/mask_search/search/           support enumeration, exact optimum, feasibility and
                                       dimension checks, and the search for W
    code/attack_complexity/            single, vectorial and distinguishing attack costs
    code/dev/                          development history, not needed for reproduction
    data/                              w24.json, w26.json (the spanning triples) and w24.log

## Links

The repository is <https://github.com/StreamCipherlinear/SNOW-Axn>.

## Requirements

* Python 3.10 or newer with numpy. The last full reproduction used Python 3.10.11 and
  numpy 2.0.2 on Windows, 2026-11-13.
* The .NET SDK or Mono for the two C# enumerations. Each of
  `code/mask_search/search/global_support_enumeration.cs` and
  `code/mask_search/search/global_support_verify.cs` defines a class with a static
  `Run()` method and no entry point, so compile it together with a one-line `Main`
  that calls `SnowSupports.Run()` or `IndependentEnumeration.Run()`.

The scripts locate `verify.py` and the data from their own path, so they run from any
directory, for example:

    python code/mask_search/correlation/verify.py

A per-script verification record is in `code/VERIFICATION.md`.

## Reproducing the results

| Paper object | Script | Expected output |
| --- | --- | --- |
| Table 1 (LAT and ILAT absolute-value distribution) and Eq. (11), the optimal and designers' triples | `code/mask_search/correlation/verify.py` | `log2 tau* = -54.7706`, designers `-55.7414`, gap `0.9709` bits |
| Section 3.3 support enumeration, Table 2 (the four supports) | `code/mask_search/search/global_support_enumeration.cs` | `eq18 all=158317600 screened=373544 feasible=4 same=4`, `lower all=314668232 screened=552044 feasible=0`, four hits |
| Section 3.3, the same enumeration by row-space membership | `code/mask_search/search/global_support_verify.cs` | the same totals and the same four hits |
| Section 3.3 byte matrices and MDS blocks (proof of the four supports) | `code/mask_search/search/global_support_matrix.py` | both matrices agree on all 256 basis pullbacks, all 276 minors nonzero |
| Section 3.3 and Table 2: exact optimum of the four supports | `code/mask_search/search/global_optimum_scan.py` | `mu=7070 ... log2tau=-54.770580312`, one optimal pair with 100 triples; three at `-55.155870468` |
| Table 3 and Eq. (12): correlation-attack cost | `code/attack_complexity/single/audit_attack_costs.py` | `bg=241 ... time=250.919215250`; designers' trail under the rule of Section 3.4, `time=252.201402061`; designers' rule at the exact trail, `time=257.186283588` |
| Auxiliary 24-mask span (development check for Section 4.2): fourth and sixth power sums | `code/attack_complexity/distinguishing/verify_list.py` | rank 24, 6,138,951 non-zero, `log2_sum4=-216.8436669675`, `log2_sum6=-327.2032268545` |
| Table 5 rows v = 25, 26 and the numbers of Sections 4.2 and 4.3: power sums of V_26 | `code/attack_complexity/vectorial/sei_v26.py` | `SEI_2(V_26) = 2^-216.797094`, `SEI_3(V_26) = 2^-327.184547`, `log2 N at v=26 = 226.0100` |
| Section 4.3: separated-time kernels | `code/attack_complexity/distinguishing/verify_mask_relations.py` | `T=(0, 2, 3) rank=383 nullity=1`, `T=(0, 1, 2, 4) rank=480 nullity=32` |
| Table 5 and the Result of Section 4.2: vectorial cost | `code/attack_complexity/vectorial/cost.py`, `bg_opt.py` | `v = 26` at the working point: `log2 N = 226.010`; `bg_opt.py`: `B_g = 214`, total data `227.857`, time `228.918`, memory `227.857`, gain `22.001` bits in time and `22.171` in data |
| Section 4.1: profiles of the four feasible supports | `code/mask_search/correlation/profile_supports.py` | optimal `11093111` pairs / `3295` levels / `-54.771`; others `11193038/3165`, `11032069/3266`, `11082491/3299`, all `-55.156` |
| Section 4.2: capacities for v = 24..29 (auxiliary) | `code/attack_complexity/vectorial/bg_scan_v.py`, `coset_v27.py`, `coset_v28.py` | auxiliary robustness scan, not quoted in the paper: `Bg=214`, data `227.857`, time `228.931/228.924/228.918/228.911/228.903/228.900`, time limit `2^228.471`, curve minimum at `Bg=214` (Appendix A, Table 7) |
| Section 4.2: supports outside the profile | `code/attack_complexity/vectorial/coset_other.py`, `coset_all.py` | each further triple adds `0.098` bits to `log2 SEI_2`; all three add `0.277` bits (data fall `0.277`) |

The remaining scripts support the argument rather than a single number. In
`code/mask_search/correlation`, `verify_exact.py` and `profile.py` enumerate the
correlation profile of Section 4.1, `levels.py` its level list, `ties.py` the
bytewise maximizers of the 24 strongest, `fourth2.py` the power sums and the check-level costing,
and `check_triple.py`, `unit_tau.py` and
`t24.py` recompute or cross-check the same numbers. In
`code/mask_search/search`, `feas.py`, `feas2.py`, `screen24.py` and
`structure_check.py` check the feasibility and dimensions of the four supports, while
`verify_space.py`, `push_k24.py`, `rank_check.py` and `crosscheck_sei.py` produce the
space W and the costs of Table 5. `pair_search.py` and `search_space.py` explore
alternative selections of the 24 triples. In `code/attack_complexity/vectorial`,
`cost.py` estimates the cost of the vectorial attack of Section 4.2 and prints the dimension scan.

## Data

`data/w24.json` holds the 24 triples of the earlier 24-mask space (auxiliary, not the W of the paper), with byte 0 leftmost; the selection
rule is recorded in the file. Its SHA-256 is
`68ef228943f9f4a2fafa571e28ce43bc06a9437b68662c12b22b8fe927a256e1`, and
`data/w24.log` is the output of the last run of
`code/attack_complexity/distinguishing/verify_list.py`; the run time is printed but not
written, so the log is byte-identical across runs.

`data/w26.json` holds the 26 triples that span V_26, the 22 of levels 0 to 7 and the four
of level 8 of the correlation profile, with the selection rule in the file. Its SHA-256
is `112361d5c62f0104bda51a64240812d47f597032b73e13fdcd6ba5a497f4bbeb`.
`code/attack_complexity/vectorial/sei_v26.py` prints the power sums and takes about four
minutes.

## Runtime

`verify_list.py` enumerates 2^24 sums and takes about two minutes, `sei_v26.py` enumerates
four cosets of 2^24 sums in about four minutes, `push_k24.py` about
one minute, and the two C# enumerations about ten seconds each including the build. The
remaining scripts finish in seconds or in tens of seconds.

## Licence

The code and data in this repository are released under the MIT licence; see
`LICENSE`.

## Citation

    @misc{snowa2026,
      title  = {Linear Cryptanalysis of {SNOW-Axn}},
      author = {{The SNOW-Axn authors}},
      year   = {2026},
      note   = {Manuscript; venue pending}
    }
