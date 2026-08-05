# Local infeasibility: complete exact-ℤ proof

Best verified: **39,022 / 39,033** (`best_agentA_39022.json`). This session reduced the
obstruction to its irreducible core and then **proved** — with exact integer solving, not
heuristics — that the entire locally-reachable solution space is infeasible.

## The result

**64 / 64 local six-selector configurations are exactly UNSAT over ℤ.** For the six selectors
that support the residual `{x_2081, x_4287, x_11368, x_13195, x_17406, x_22562}` (all other 250
selectors held at agentA), every one of the 64 on/off patterns was decided:
- 27 UNSAT already at the linear level (the p-fractional slack obstruction), and
- 37 UNSAT via exact nonlinear solving: z3 solves the linear integer subsystem → integer
  particular solution `x0`; `flint fmpz_mat.nullspace` gives the integer kernel `K`; substitute
  `x = x0 + K·t` (integrality automatic) and z3 decides the residual degree-2 system in the
  53–88 free parameters `t`. Every combo decided in 5–51 s. 0 SAT, 0 timeout.

## Why every family of techniques fails (each checked this session)

- **Continuous repair (Newton / tangent / null-space):** the residual is a single conserved
  functional (obstruction dim exactly 1), inconsistent in free-input space, core-augmented
  space, and bit space. Refuted.
- **Non-degenerate core:** provably decoupled from the verifier (perturbing core knobs moves
  0/16 leaf-ripple residuals). Refuted.
- **Algebraic (XL / cube / HOD / elimination):** the 256-bit message codeword is dense,
  degree ≥ 3 in the bits, with cryptographically random constants; bit-Jacobian rank 256
  (message uniquely pinned); degree-2 XL is under-determined (19013 eqs vs 32897 monomials,
  ratio 0.578). No linearization shortcut.
- **Small moduli (mod-2 / 2-adic Hensel / mod-q):** the obstruction is a mod-`p` divisibility
  of a sub-`p` gap absorbed only by a `p`-granular slack; since `p` is odd and coprime to
  `q=6672769`, the slack spans everything mod 2ᵏ or `q`, so those projections are **structurally
  blind** to it. Mod-2 circuit-SAT: selector backbone empty (0/256 forced). Ruled out.
- **Exact-ℤ (HNF + z3):** agentA branch UNSAT (even with expanded compensator pool); branch B
  (1,1) linear-feasible with 20 free params but its 42-equation gadget UNSAT; all 64 local
  combos UNSAT (above).

## What remains

Crossing the final 11 equations requires a **non-local** reconfiguration of the 256-bit
message — the setter's global witness — which is a high-density (density ≈ 1) knapsack / CVP
over the random pin constants, coupled at degree ≥ 3. That is the trapdoor; no polynomial
attack on it was found, and none is expected for random constants at density 1.

## Provenance
Four independent agents (cube/HOD, core-steering, algebraic-elimination, modular/SAT) plus the
main analysis converged on this, with complementary proofs. Scripts: `local_tangent.py`,
`closure_test.py`, `conserved*.py`, `bit_newton.py`, `reaudit_*.py`, `z_*.py`, and the SAT
agent's scratchpad (`decide20.py`, `reduceB2.py`, `scan*.py`).
