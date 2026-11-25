# reduced_scale — reduced-scale analogue of SNOW-Axn for piling-up validation

Self-contained scripts (numpy only, no data files) for the FSE reduced-scale
experiment: an 8-bit analogue of SNOW-Axn whose Monte Carlo / exact
correlations are compared with the piling-up predictions of the vectorial
attack (linear approximations, folded checks, capacity).

Requirement: **Python 3.10** with numpy (the paper scripts use `python3.10`;
the default `python` on this machine is AutoClaw's 3.13 without numpy).

## The analogue

* words are bytes, FSM registers `R1,R2,R3 ∈ GF(2)^8`, `S` = AES S-box,
  `σ = rotl3`, `L = I` (so `E = S`);
* LFSR: eight bytes `u_0..u_7` with the design's byte recurrence
  `u'_r = c_r·u_r ⊕ u_{r+1} ⊕ u_{r+3}` (indices mod 8), `c_r = α` for `r<4`
  and `β` for `r≥4`, where `α,β` are roots of `x^8+x^5+x^4+x^3+1` and
  `x^8+x^7+x^5+x^2+1`;
* `H(t)=u_0(t)`, `z(t)=R1(t)⊕R2(t)⊕H(t)`,
  `R1(t+1)=σ(R2(t)⊕R3(t)⊕H(t+1))`, `R2(t+1)=S(R1(t))`, `R3(t+1)=S(R2(t))`;
* three-clock linear approximation with masks `a,b,c` (cf. Eq. (N) with L=I):

```
e(t) = a·[z(t-1)⊕H(t-1)] ⊕ b·[z(t)⊕H(t)] ⊕ c·[z(t+1)⊕σH(t+1)⊕H(t+1)]
     = a·[S^-1(R2)⊕S^-1(R3)] ⊕ b·[R1⊕R2] ⊕ c·[σ(R2⊕R3)⊕S(R1)]
```

and, with `c' = σ^-1(c)`, the **verified** piling-up correlation

```
tau(a,b,c) = LAT[b][c] · LAT[a][b^c'] · LAT[a][c'] / 2^24 ,
```

where `LAT[a][b] = sum_x (-1)^{a·x ^ b·S(x)}`.  Note: this differs from the
expression in the work order, which prints `LAT[b][c']` as the first factor.
`verify.py` V8 checks by exhaustive enumeration of the 2^24 register states
that `E[(-1)^{e(0)}] = tau` for the first form and that the printed variant
does not.  For the mask convention of Eq. (N) (`C` applied to
`sigma^-1(z^{t+1} ^ sigma H3 ^ H2)`) the verified expression is
`LAT[b][σ(c)]·LAT[a][b^c]·LAT[a][c]/2^24` (verify.py V8b).

## Files

| file | purpose |
|---|---|
| `build_tables.py` | GF(2^8) / AES S-box (computed from scratch) / α, β / σ / LAT tables |
| `verify.py` | structural checks V1–V9 (identities, LFSR/FSM rank/bijectivity, two forms of e(t), exhaustive t=0 piling-up, both mask conventions) |
| `search_masks.py` | exhaustive search of all `255^3` mask triples, top list, independent basis, span; writes `masks.json` |
| `monte_carlo.py` | Monte Carlo over independent initial states (N = 2^22 … 2^26), batch-mean standard errors |
| `exact_folded.py` | exact (zero-MC-error) folded correlations and capacity by a 2^24-state recursion, with brute-force checks |
| `collect_results.py` | reads the JSON outputs and prints the report tables |

## Run

```bash
python3.10 verify.py                                   # ~4 s, must print ALL CHECKS PASSED
python3.10 search_masks.py --top 40                    # ~1 s
python3.10 monte_carlo.py --n 24 --batches 32 --json mc_n24.json    # ~13 s
python3.10 monte_carlo.py --n 26 --batches 32 --json mc_n26.json    # ~47 s
python3.10 exact_folded.py --check --json exact_results.json        # ~6 s
python3.10 collect_results.py                          # Markdown tables
```

Outputs:
`mc_n24.json` / `mc_n26.json` (single-clock, folded d=1,2,3, independence test,
capacity), `exact_results.json` (exact folded + capacity), `masks.json`
(masks/span), and the human-readable tables from `collect_results.py`
(same columns as `_rev/fse2/report_tables.md`).

## Fixed-format output conventions

* `search_masks.py`: header comments then one row per triple
  `rank a b c c' tau_num tau log2|tau|`.
* `monte_carlo.py`: sections `(i)` single clock, `(ii)` folded, `(iii)`
  different-clock independence, `(iv)` capacity, then a summary; all means are
  batch means and every row carries a standard error.
* `exact_folded.py`: one row per span element with
  `tau_pred tau_exact C2(1) C2(2) C2(>=3) C2/tau^2 bits`, then the capacity.

## Notes / caveats

* The toy LFSR transition matrix has rank **61** (not 64; a 3-dimensional
  kernel), so H values are linearly dependent over long spans
  (rank(H(0..7)) = 37).  Every H-window of length ≤ 4 used in this experiment
  is nevertheless full rank (`8d` for d ≤ 4), which is all the exact folded
  computation needs.  This is a property of the *fixed toy parameters*, not of
  the real 512-bit LFSR.
* The folded correlation of two clocks separated by `d` is
  `C2(d) = E[(-1)^{e(t) ^ e(t+d)}]`.  In the analogue it equals `tau^2`
  exactly for `d = 3` (and `d = 4`), i.e. when the two three-clock windows are
  disjoint, and deviates for `d = 1,2` (overlapping windows).  For `d = 1` the
  capacity is 9.73 bits above `Σ tau^4`.  Consequence for the paper: state the
  folded offset, and prefer `d >= 3`.  Offsets `d >= 5` are not covered: the
  toy LFSR has linearly dependent H bytes over five-clock spans.
* The Monte-Carlo capacity estimator is noise-dominated at N = 2^24/2^26
  (the bias of squaring a correlation estimated to 1/sqrt(N) is much larger
  than Σ tau^4 = 2^-34.4); the capacity is therefore reported from the exact
  state-space computation, with the individual folded correlations validated
  against the Monte Carlo (chi^2/ndf = 0.91 over 21 values).
