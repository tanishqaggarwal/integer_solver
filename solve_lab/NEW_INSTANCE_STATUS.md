# New instance (`EQUATIONS.txt` @ b616370) — progress & obstruction

Honest status: the re-randomized instance is the **same trapdoor family** but its twist is
**materially harder** than the previously-solved one. I reconstructed the partial
generically and isolated the twist, but have **not** cracked the final coupled system.

## What works (verified)
- Regenerated atoms (`poly_atoms.py`) and gates (`extract_atoms.py`): same family —
  degree histogram `{1:20124, 2:25480, 4:694}`, 841 huge-constant load atoms, a 220-var
  identity wire, ~694 perfect-square verifier atoms.
- **`rebuild_partial.py`** — a *generic* reconstruction (greedy topological gate
  orientation + forward-eval from free-inputs = 0, the same principle proven on the old
  instance where every free input is 0). Reaches **39,005 / 39,033** and isolates the twist
  to **4 nonzero atoms**. This confirms the core methodology transfers.

## The twist (4 atoms) and why it is harder
Two product-slack chains, exactly analogous to the old F/H:
- Chain 1: `x_14257 = x_7497·x_23917` must absorb `8863713·(x_18956 − BIGCONST)`.
- Chain 2: `x_32989 = x_11436·x_22399` must absorb the `HUGE2` gap.
- x_7497, x_11436 (and deeper x_22820, x_14393) are genuinely free (1 atom each).

**The decisive difference from the solved instance:** there the 220-var wire was a *free*
parameter, so I set it to 1 and the gaps absorbed as `1 · (−G)` — a trivial factorization.
Here the wire multipliers x_23917, x_22399 are **pinned to a fixed V0 ≈ 2²⁵⁶** by a hard
constraint (`x_26064 − V0` is a real shared atom, and a 5-hop identity chain ties
x_22399/x_23917 to x_26064). With the wire pinned:
- the product slacks can only produce **multiples of V0**;
- the gaps (BIGCONST, HUGE2) are **not** multiples of V0 (verified);
- so the non-V0 part of each gap must route into the two shared verifier checks
  **40907, 44255**, where the two chains' huge values must **cancel**.

The cancellation is not linear: `a·BIGCONST + b·HUGE2 (mod V0)` is ~2²⁵⁶ in both checks
(not reachable by the small linear vars), and the checks contain **products of
already-huge variables** (e.g. `x_29356·x_33469` with x_33469 ≈ 2²⁵⁶). Absorbing the
residue therefore requires solving a **nonlinear (product) Diophantine system** over ℤ —
the pinned wire has turned the old unit-factorization into a residue/factoring-type
problem. This is plausibly the intended hard core of this instance.

## Approaches tried (and why they stalled)
- `local_solve.py` — Smith-normal-form integer solve over spine+checks. Solves the spine,
  but the ℚ/linear treatment either blows up nonlinear check terms or ripples into
  neighboring atoms (39,005 → 38,963/38,992).
- `auto_solve.py` — auto-expanding frontier + integer solve. The coupled component grows
  without bound (the checks couple to thousands of vars), so a single SNF does not
  converge in memory/time.

## Honest assessment / next angles
The remaining core is a small but **nonlinear, V0-modular** integer system. Cracking it
likely needs one of: (a) the intended gate orientation that keeps the huge flow off the
small-coefficient checks (if one exists); (b) a targeted nonlinear solve that chooses the
product multipliers (x_29356, …) to hit the two check residues mod V0 while the V0-quotient
is handled by the free partners; or (c) recognition that this instance encodes a
factoring-hard verifier, in which case it is not solvable by local repair without the
trapdoor secret. Progress and tooling are committed; the 39,005/39,033 partial is real.
