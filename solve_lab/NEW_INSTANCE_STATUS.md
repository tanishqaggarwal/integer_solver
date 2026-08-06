# New instance (`EQUATIONS.txt`) — progress & the reduced core

Honest status: **best verified = 39,013 / 39,033** (`best/new_instance_partial_39013.json`).
The 20 remaining are the quadratic verifier squares — the trapdoor core. This session
**fully reduced that core to three scalar conditions** and localized them to two deep
control gates, but the final integer solve is blocked by global wiring coupling.

## What works (verified)
- `forward_construct.py` — topological forward-construction with double-width split. Solves
  ALL linear wiring (39,013/39,033). The only failures are the 20 quadratic squares.
- `checker.py best/new_instance_partial_39013.json` → 39,013/39,033. Solid.
- Scale: 8,583 free inputs, 30,165 gate outputs; ~10,517 genuinely-constraining equations
  (28,516 are auto-satisfied definitions). Quadrant (1,1), activators (24601, 2081).

## THE REDUCTION (this session's breakthrough) — see CORE_REDUCTION.md
All 20 core equations are exact integer combinations of three "monster" quantities:
- **M1** = x_15298·x_11150 + x_4007          (load L1 = x_11150)
- **M2** = x_15298·x_25739 − 6672769·x_29804  (load L2 = x_25739)
- **M3** = 537773·(x_15298·x_37758) − x_35605 (load L3 = x_37758)

Core ⟺ **M1=M2=M3=0**. The three loads are wiring-defined by two base gates
S=x_35389, T=x_6671 (coeff matrix has invertible 2×2 minors), so:

  M1=M2=M3=0  ⟺  **S ≡ 0 mod p  AND  T ≡ 0 mod p**  (+ one M2 mod-6672769 cond)

and S,T collapse to two control differences (each with a free-input side):
- **x_29322 = x_14853 − x_12186**  (both free)
- **x_3558  = x_24908 − x_16742**   (x_16742 free)

with S ≡ −x_3558² and T ≡ −x_3558·x_1326  mod p (given x_29322≡0). So the ENTIRE core is:
**x_3558 ≡ 0 mod p and x_29322 ≡ 0 mod p**, then set private quotient handles
x_30317, x_2936, x_5146 for M1, M3, M2.

## VERIFIED partial results
- Minimal residue fix (x_14853 −= x_29322%p, x_16742 += x_3558%p, both < p) drives
  **S ≡ T ≡ 0 mod p** and clears **M1, M3 exactly** (L1%p = L3%p = 0). Confirmed.
- The fix breaks ~23 non-core wiring equations. Their Jacobian is rank-9 and CONSISTENT
  (solvable), but the spanning handles are high-coupling: every heal attempt (greedy Dixon,
  constrained residue-pinned, ripple-minimized) diverges via ripple. Coupling closure is
  ~the whole system (21,825 eqs / 5,878 handles).
- M2 residual: L2/p % 6672769 = 645924 ≠ 0; tunable via the quotient parts (d,m) of
  x_29322, x_3558 (one quadratic condition mod 6672769, 2 knobs → solvable in principle).

## The remaining obstacle
Zeroing x_3558, x_29322 residues requires moving within the wiring's null space. The
control vars are multi-role (x_14853 also feeds x_9192; x_16742 also feeds x_27713), so
naive residue fixes ripple globally. Need either (a) a proper global mod-p null-space move,
or (b) a from-scratch construction that bakes S≡T≡0 in without losing handles for the ~26
equations that currently consume x_14853/x_16742.

## Tools added this session
diag_terms, diag_monsters, diag_coupling, diag_cone, diag_ctrl, diag_consumers, diag_rank,
diag_private, diag_closure, diag_constraining; exp_minfix, exp_quot, exp_footprint;
fc_zero, heal_dixon, heal2..5 (all diverge — documented). CORE_REDUCTION.md is the summary.
