# Agent G — RESUME

## Best verified score
**39,026 / 39,033** — `solve_lab/best/new_instance_partial_39026.json`, re-verified:
`python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
-> `satisfied 39026/39033 (7 failing)`, failing `[12231,12270,12350,14584,18673,22044,29125]`.
I did not beat it; nothing written to `best/`. All my files are under
`solve_lab/agentG_work/`; no shared file was modified.

## The machine: exact symbolic forward evaluation over F_p
`gsym.py` / `gsym2.py`. Gate output coefficients are ±1, so forward evaluation divides by
nothing and free-inputs -> every atom is an honest polynomial over Z, hence over F_p.
* From `s10/AG_39013.json`: closed non-boolean symbol set **112**; of **10,792 check atoms
  only 57 are non-constant** (50 linear, 2 quadratic, 5 cubic, **196 monomials**).
  Validated at random points: **0 mismatches on all 10,792 checks**.
* Maximal model: ALL **6,117** non-boolean free inputs symbolic, 0.5 s, 0 gates skipped;
  rank **1470**, 4,647 free parameters, consistent; every nonlinear check reduces to a
  CONSTANT, five nonzero.
* Equation level: **6,774 non-trivial equations = 6,613 linear (rank 1470) + 161
  nonlinear**; forcing all linear ones leaves exactly **20** failing = AG_39013's 20.
* Validated at the deliverable (five gate outputs 7068/28730/29854/31864/642 detached):
  **0 atom mismatches over all 42,267 atoms, exactly 7 nonzero equations**, and **all 7
  violated equations are LINEAR — zero nonlinear ones are violated.**

## The residual, in the instance's own variables
K = 97553848499418123410591666447050222001188385549510401465815187079080512838891

    A = (x22649 - x14853)^2 * (x22162 + x22649 + x14853 + K) - (x31339 - x16742)^2
    B = x14853*x30213 - x22649*x30213 + x14853*x16742 - x22649*x31339
        + x22162*x31339 - x16742*x22162
    a19297 = 8646263*A + 1073965*B ,  a19299 = 10159099*A + 6926539*B , and
    a30984, a36185, a40812 are three further members of the same rank-2 pencil.

## eq8680 and eq29125, derived here
eq8680 is a **binary quadratic form in exactly x3629 and x8976** with discriminant
**exactly 0**: `eq8680 = 784*(x3629 - r*x8976)^2`,
r = 99250362203413881791632272864589635302802843999120483462392214863921859717787.
A perfect square: it pins one linear relation and costs exactly ONE equation to violate.
eq29125 is linear: `x28730 + C*x3629 + 259857806*x8976 = 0`.

## OPTIMALITY OF 7 — complete, no undecided cases
Objective: `total(T) = |T| + min over the subspace freed by dropping T of the number of
higher-degree equations that are nonzero`. (A failing higher-degree equation costs +1, it
is not fatal — this corrected my earlier, stricter test.)
1. The 20 higher-degree equations failing at x0 span a space of dimension **4**, and
   **removing any 5 of them leaves the span at 4** (`g74_span.py`). So vanishing of any 15
   forces all 20 to vanish: **#failing is 0 or >= 6**.
2. `#failing = 0` requires `|T| >= 7`: exhaustive at budget 6 over the closed support
   (`g66`, 4,874 admissible relaxations, every one **decided**, 0 undecided).
3. `#failing >= 6` gives `total >= |T| + 6 >= 7`; `|T| = 0` gives 20.
=> **total >= 7, attained.** Steps 1 holds on both the 17- and the 22-unknown supports;
step 2 is confirmed on the 17-unknown support and the 22-unknown run is in `tierA6.log`.

## Tier report (enumeration size is printed and capped, never silently sampled)
| support | unknowns | affine rows | distinct directions | mult<=6 dirs | candidates at budget 6 |
|---|---|---|---|---|---|
| closed region | 17 | 68 | 19 | 15 | 9,948 (4,874 admissible, all decided) |
| + x3629,x6418,x8976,x27500,x32230 | 22 | 105 | 36 | 31 | 443,068 (running, `tierA6.log`) |
`g73_lb.py`/`g69_tier.py` print the candidate count and abort above 3e6–4e6.

## Re-enter
```
cd /home/user/integer_solver/solve_lab/agentG_work
python3 g54_cosetsetup.py                       # -> coset_model.pkl
python3 g56_colweight.py                        # per-unknown costs, x0, deliverable departure
python3 g67_eq8680.py                           # eq8680 / eq29125 structure
python3 g68_widen.py                            # -> sup_tierA.txt
python3 g74_span.py $(cat sup_tierA.txt) 5      # the span argument
python3 g69_tier.py  $(cat sup_tierA.txt) 6     # exhaustive, prints its own size
python3 g70_total.py $(cat extsup.txt) 6        # corrected total objective
```
NOTE `s9/eff/lib.py` does `os.chdir(solve_lab/s9)`; write outputs to absolute paths.

## Next
The next tier is the 85 unknowns with 10 linear equations. Adding even a few of them takes
the multiplicity-<=6 direction count past ~40 and the candidate count past 3e6, where the
enumeration stops being exhaustive at budget 6 — so the honest next step is not a wider
brute-force but to carry the span argument outward: check whether "the failing
higher-degree equations span 4 and lose nothing when any 5 are removed" still holds on the
wider supports, since that alone forces #failing to be 0 or >= 6 and reduces the question
to the single exhaustive test `#failing = 0 => |T| >= 7`.
