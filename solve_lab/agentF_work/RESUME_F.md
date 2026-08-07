# RESUME_F — agent F (multi-modular / p-adic lifting).  NO curve/group framing; integer polynomials only.

## READ THIS FIRST — the infeasibility argument is WITHDRAWN.  The instance is OPEN.
One thing is proved: **rank(M) = 39,033, dim ker(M) = 0** (three independent computations).  That closes
the "cancelling nonzero residual" route and makes all-atoms-zero *equivalent* to a full solve.
It does NOT make the instance infeasible.  My section-7 exhaustion assumed at most one ON boolean per
OR-tree ("link (i)").  A rigorous mod-p rigidity engine (`cfg_rigid2.py`) now shows that assumption is
UNSUPPORTED: on same-tree pairs it derives no contradiction (45 pairs, 0 conflicts), both booleans' pin
wires stay forced to their own constants, and the wire that lost its forcing is an ADDER
`x11317 := x11532 + x14681` that carried the whole selected value in the single-boolean case.  That is an
accumulator signature.  A numerical run with two same-tree booleans also activates a SECOND gate
(x24533) with three further checks.  So the reachable set of the selected pair may be far larger than the
13,884 one-per-tree combinations I enumerated.  **No infeasibility is claimed by me.**

## PROVED:  rank(M) = 39,033 and dim ker(M) = 0  (exact, certificate-verified)
M = the 39,033 x 39,033 equation-atom incidence matrix (525,982 nonzeros, coefficients from the spine
decomposition).  A characteristic-free peeling cascade -- a row whose surviving support is a single atom
forces that atom to 0 -- starts from the unique degree-1 row and **forces all 39,033 atoms to zero**.
The elimination order is stored as a checkable certificate (`peel_order.npy`) and re-verified by an
independent pass (`peel_cert.py` -> `certificate verified: True`).  Holds over Z and over every field of
characteristic > 80.
Independently confirmed by **Wiedemann over two word primes**: q = 2^31-1 gives minpoly degree 39,033 = n
with trailing coefficient 268435456 != 0, and q = 2147483629 gives degree 39,033 with trailing coefficient
11716781 != 0 -- M nonsingular both times (implementation validated first on 4 nonsingular and 3
rank-deficient controls with independently computed true ranks).  The pivots of the peeling cascade are all +-1 or +-2, none divisible by
any odd prime tested, and three randomized peeling orders all force 39,033/39,033 atoms.
**Consequence: any assignment satisfying all 39,033 equations must make all 39,033 atoms exactly zero.**
The all-atoms-zero model is therefore not a restriction but an equivalence; the "cancelling nonzero
residual" route to a full solution does not exist, and every frame-optimality result in this lab that was
conditional on that model loses that condition.

## Scores
- Shared baseline **39,026** re-verified by me (`solve_lab/best/new_instance_partial_39026.json`).
- My own pipeline's best: **39,024** = `agentF_work/best_F_39024.json` (checker-verified, 9 failing).
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
1. **Settle the accumulator question** -- this is now the top item and it is a CONSTRUCTIVE opportunity,
   not a barrier.  With two same-tree booleans ON: does the selected wire equal the SUM of the two
   per-boolean contributions?  If yes, the selected pair ranges over subset sums of the 256 pin constants
   and the enumeration that produced "no solution" was over a tiny slice of the real space.
   Tools ready: `cfg_rigid2.build(bits)` (~15 s, exact, 0 conflicts everywhere so far) and `gs2.solve`.
   The right experiment: take two same-tree booleans, and instead of letting `gs2` zero the contributions,
   force the selected wire to C_j1 + C_j2 (mod p) and see whether the residual atoms close.
   Also characterise the extra gate x24533 that turns on with two same-tree booleans and its 3 new checks.
2. Finish `sweep_ii.py` (exhaustive link (ii), ~15 s/boolean, checkpointed to `sweep_ii.json`).
3. Continue `sweep_i.py` only as a diagnostic -- it is currently inconclusive by construction, so treat
   its coverage percentage as "fraction tested", never as "fraction excluded".
2. Coset search for a better cancellation frame: the number of the 12 frame rows that cancel is decided by
   the coset of the residual value in the 7-generator lattice (invariants: r1 mod 7376877, r2 mod p,
   r3+r4 mod p, r5-r6 mod p).  Sweep the ON-boolean pair and the break placement to move the coset; a
   6-row coset would give 39,027.  (Note: r2 mod p can never be 0, since the two boolean groups' pin
   constants are disjoint mod p, so equation 29125 appears to be unavoidable in this frame family.)
3. mod p^k with a proper Hensel step rather than greedy repair, to characterise the mod-p residue variety.

## Multi-modular results (new)
- FULL system solved (0 broken atoms, 0 failing equations) modulo: all primes < 110, 1009, 10007, 100003,
  1000003, 10000019, 10^9+7, 2^31-1, 2^61-1, 2^89-1, 2^127-1, 2^255-19; prime powers 2^64, 2^100, 3^40,
  3^80, 5^10..5^40, 7^25, 11^20, 13^20, 1009^8, 65537^4, 1000003^3, (2^31-1)^2; and the handle multipliers
  M themselves.  Checkpoints in `modm_results/`.
- CRT of 20 independent 7-digit-prime solutions gives a solution modulo a **399-bit** composite Q with
  0 failing equations mod Q (`crt.py`, vectors in `crt_sols/`, report `crt_report.json`).
  Its balanced integer representative scores 38,991 (`crt_balanced.json`) -- CRT alone does not lift.
- Handle structure validated at scale: 3,173 genuine divisibility atoms all have handle ≡ 0 mod p under
  4 random draws of all free inputs.
- **Obstruction is at exactly one modulus, p, at level p^1 only.**
- **mod p and mod p^2 do NOT solve**: same code, same booleans -> 21 nonzero atoms / 124 failing eqs (p) and
  15 / 102 (p^2), versus 0/0 in 3.5 s for all 60 other moduli.  Checkpoints `modm_results/P_p.json`,
  `modm_results/P_pp2.json`.
