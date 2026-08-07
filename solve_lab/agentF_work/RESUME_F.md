# RESUME_F — agent F (multi-modular / p-adic lifting).  NO curve/group framing; integer polynomials only.

## Scores
- Shared baseline **39,026** re-verified by me (`solve_lab/best/new_instance_partial_39026.json`).
- My own pipeline's best: **39,024** = `agentF_work/F_best.json` (checker-verified, 9 failing).
  Also `best_F_39022.json` (39,022, only 2 nonzero atoms) and `F_frame.json` (39,023).
- Nothing above 39,026 yet.

## Established (mine, independent)
- 39,033 atoms; 30,001 definitions form a **DAG**; 8,747 free inputs; **9,032 residual atoms**.
  **All residual atoms zero => all 39,033 equations hold.**
- All-free-inputs-zero => **39,005**, only 3 nonzero atoms.  Engine `fwd.py`, 16 ms/pass, matches checker.
- Every check has the exact form `A - B = M*p*h` with h a free integer, M a ~24-bit literal, and
  p = 115792089237316195423570985008687907853269984665640564039457584007908834671663 a 256-bit prime
  literal of the file (it defines x26064).  So every check is a congruence mod M*p.
- **The lift obstruction concentrates at the single modulus p.**  `modq.py` fully solves the entire
  39,033-equation system mod q (0 broken atoms, 0 failing equations) for 13 primes q != p including
  2^127-1 and 2^255-19 (size control).  Reason: `M*p*h` is surjective mod any q coprime to M*p.
- **No obstruction above p^1**: with the two unconditional constant pins relaxed, everything closes exactly
  over Z (`relaxed_pin.json`) -- handles absorb every quotient.  A mod-p solution would lift immediately.
- All-atoms-zero model is exhausted (13,884 combinations, zero solutions) but **no infeasibility is claimed**:
  a solution may carry nonzero atoms that cancel; that needs ker(M) != 0 for the 39,033x39,033 incidence
  matrix M (525,982 nnz), which I have NOT computed.
- 39,026 optimality in its own frame independently CONFIRMED (unique integrally-reachable 5-subset of 12
  rows; all 924 6-subsets solvable over Q -- integrality is the barrier).

## Re-enter
    cd /home/user/integer_solver/solve_lab/agentF_work
    python3 parse3.py; python3 circ4.py; python3 sched.py; python3 supp.py   # rebuild pickles (~1 min)
    python3 fwd.py       # 39,005 from all zeros, 3 broken atoms
    python3 modq.py      # full solve mod q for many q
    python3 frame.py     # frame/lattice analysis helpers ; intsolve.py = exact integer HNF solver

## Next experiments (in order)
1. Prime POWERS: solve mod q^k for q != p and mod p^k; confirm Hensel lifting is free away from p.
2. CRT-reconstruct a solution mod a large composite Q = prod q_i and measure how far it is from an integer
   solution (balanced representatives) -- quantify the size gap.
3. Compute rank(M) (or at least test ker(M) != 0) -- that is the only remaining gate on the all-atoms-zero
   argument, and the only route to a residual that cancels instead of failing.
