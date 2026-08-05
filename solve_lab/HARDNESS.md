# Fundamental hardness after local infeasibility (working note)

Channel (1,0) is proven ℤ-infeasible (see INFEASIBLE_CHANNEL_10.md). The residual
hardness is finding the setter's NON-LOCAL selector witness.

## Search structure (measured)
- Master select `x_15298 = OR(tree1: 178 bits) · OR(tree2: 78 bits)` — 256 free
  boolean selector bits. best_agentA has Hamming weight 2 there: {x_2081(=a, in tree2),
  x_24601}. 844 boolean free inputs total; 8583 free inputs total.
- mod-p rank: 3036 active inputs, 5548 dead (mod-p-inert, multiply the hardcoded p-wire).

## Why the cheap square-root route is DEAD
- The message load is additive (subset-sum over 41 bits) — the classic MITM target —
  BUT all 41 load bits are **mod-p-inert for the obstruction**: flipping any leaves
  G1,G2 unchanged mod p (verified, hardness_loadbits_inert.py), only breaking other
  positions. No 2^(n/2) knapsack shortcut on the obstruction.
- Single-bit flips of the 254 non-active selector bits are also mod-p-inert for Q
  (hardness_selector_flat.py); only local a=x_2081 moves Q, and its 4 combos are UNSAT.
- ⇒ Q is **flat/avalanche to single-bit changes** — a preimage signature, not a soft
  additive target. The obstruction is dim-1, reciprocal-locked (r1·r2 ≡ 1 mod p).

## Open: exact exponent + square-root opportunity (under analysis by 2 agents)
- Naive: search the ~256-bit selector witness for a mod-p-feasible (and liftable) config.
  Not the inert subset-sum; a preimage-like search ⇒ ~2^128-class unless structure helps.
- Genuine sqrt candidate: the reciprocal lock r1·r2≡1 and the two-tree PRODUCT
  OR(178)·OR(78). IF the joint solvability condition factors as A(left bits)·B(right bits)≡1
  over disjoint halves → birthday collision ~2^(n/2). Decisive question being tested:
  does r1 / the message-vs-response condition split across disjoint variable halves?
