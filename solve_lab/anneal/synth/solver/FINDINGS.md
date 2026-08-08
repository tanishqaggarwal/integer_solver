# Does the encoding anneal? — and what that does to the run count

All measured on synthetic planted-key instances (`synth/gen.py`), classical Ising
solvers standing in for the QPU. Code: `model.py` (QUBO→numpy Ising), `solvers.py`
(SA, parallel tempering, tabu, simulated bifurcation), `diagnose.py`, `bench.py`,
`restart_scaling.py`.

## 1. No solver beats plain SA, and all stall early

One modular multiplication `a·b ≡ c (mod p)`, free search, mid budget, 5 reps —
`best E` is the lowest energy reached (0 = solved):

```
  s    n     SA        PT        tabu      SB
  6   153   0/5 E1    1/5 E0    0/5 E1    0/5 E8
  8   255   0/5 E4    0/5 E2    0/5 E1    0/5 E18
 10   359   0/5 E5    0/5 E3    0/5 E9    0/5 E23
```

Parallel tempering and simulated bifurcation — the two algorithms that beat SA on
*barrier*-dominated problems — do **not** beat it here. Past ~6-bit words nothing
reaches the ground state at any budget tried. That already rules out "there is a
better annealing schedule."

## 2. The landscape has no gradient — the decisive measurement

`diagnose.py` fixes every ancilla to its forced value and measures energy vs
**digit-distance to the unique solution** (so E reflects only the verifier's signal):

```
bits=16 mu=10:  dist   minE  meanE  maxE
                  0      0     0.0    0        <- the solution
                  1     15    19.2   22
                  2     12    17.6   24
                  3     11    18.0   26
                  4     10    18.2   27
                  5     11    18.0   26
   correlation(distance, energy) = 0.060
```

A candidate that differs from the key in a **single digit** has essentially the
same energy (15–22) as one that differs in **every** digit (10–27). There is no
slope toward the solution: it is an isolated hole (E=0) in a flat plateau. This is
a *needle in a haystack*, not a rugged funnel — confirmed by §1, where the
barrier-crossing solvers gained nothing.

`restart_scaling.py`: even at μ=2 the SA hit-rate per restart is <0.7% (1/p > 150),
already far worse than the 2^μ=4 a pure needle search would give — because the
arithmetic ancillas also do not anneal, stacking on top of the flat answer space.

## 3. What this does to the number of runs — the honest bottom line

The interval split gives **outer runs = 2^(b−μ)**, and shrinking the encoding
raises μ, cutting outer runs. That part is real. But with no gradient a single
anneal finds a sub-instance's needle only by chance, so **anneals per sub-instance
≈ 2^μ** (empirically worse). Therefore

> **total anneals = 2^(b−μ) · 2^μ = 2^b, independent of μ.**

Raising μ (shrinking the window) trades outer runs for inner difficulty at a fixed
product. The annealer delivers **no speedup over classical brute force**, and none
of the encoding minimisation — real as it is, ~3× per window, at the arithmetic
floor — changes that. It changes only whether one sub-instance *fits*; it cannot
make the sub-instance *solvable* in fewer than ~2^μ tries.

This sits on top of the size wall (`squeeze/FINDINGS.md`: one 256-bit modmul needs
≥4× a 4,400-qubit machine in partial-product ancillas alone). Even on a machine
large and precise enough to hold a sub-instance, the run count is the classical
search cost. For a b-bit synthetic key small enough to fit one QUBO on
full-precision hardware (b ≲ 44), recovery takes ~2^b anneals — the same order as
enumerating the key classically.

**Reducing the number of runs below the classical search bound is not possible for
this problem on an annealer, because the encoded landscape carries no signal.**
The clean negative is the result.
