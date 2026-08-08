# The arithmetic-annealing frontier F (measured)

F = the largest field size s (bits) whose modular multiply `a·b ≡ c (mod p)` a real
Ising solver drives to the ground state (E=0). Best solvers only (parallel tempering
`pt`, tabu); SA and simulated bifurcation are strictly worse.

## Success rate vs field size (free search, best solvers)

| s (bits) | vars | PT hits | tabu hits | reliably solved? |
|---|---|---|---|---|
| 4 | 73 | 40/40 (100%) | 40/40 (100%) | **yes** |
| 5 | 107 | 2/20 (10%) | 10/20 (50%) | marginal |
| 6 | 153 | 8/40 (20%) | 2/40 (5%) | marginal |
| 7 | 231 | 0/20 | 0/20 (best E=2) | **no** |
| 8 | 255 | 0/20 | 0/20 (best E=1) | **no** |
| 10 | 359 | 0 | 0 (best E≥3) | no |

**F ≈ 5–6 bits.** Reliable at s≤4, a coin-flip at s=5–6, dead by s=7 — no solver
reaches E=0 at any budget tried. Parallel tempering and simulated bifurcation (the
barrier-crossing algorithms) do **not** extend F: at s≥7 they stall at the same
nonzero energies as SA, so the wall is absence-of-gradient, not barrier height.

## Landscape at the frontier (energy vs distance-to-solution, modmul)

correlation(distance, energy): 0.14 (s=5) → 0.11 (s=8) → **−0.03 (s=10)**. A
one-bit-wrong operand has mean energy indistinguishable from a maximally-wrong one;
the signal is gone by s≈10. The faint signal at s≤6 is why tabu occasionally wins there.

## Consequence

The modular multiply is itself gradient-free past ~6 bits — the same needle structure
as the full DLP, and for the same reason (multiplication is a pseudorandom function of
its input bits). Every comb window does full **field-size** arithmetic regardless of how
few scalar bits it resolves, so no annealer call on a field larger than ~6 bits ever
succeeds. **F ≈ 6 is nowhere near the 32 a 32-bit key would need**, and it cannot be
pushed there: the obstruction is structural, not effort.
