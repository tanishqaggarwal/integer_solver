# RESUME_E — agent E checkpoint (session stopped by coordinator, CPU pressure)

**WITHDRAWN ANGLE — DO NOT RESUME:** generator/authorship reverse-engineering (PRNG hunting,
emission-order forensics, "recover the generator and run it forward").  Nothing below depends
on it; everything is computed from the equations as mathematical objects.

## Verified scores
- `../best/new_instance_partial_39026.json` = **39026/39033**, fails
  [12231,12270,12350,14584,18673,22044,29125].  CONFIRMED with `../checker.py`.
- My own best state: **39,015** (18 failing eqs).  Never beat the baseline.  No file of mine
  is >= 39,026, so none was promoted.

## Rebuild the pipeline from source (in this directory, ~1 min total)
    python3 parse3.py          # -> model3.pkl   (40,727 atoms; 9,710 eqs are squares S*S)
    python3 dag.py             # -> dag.pkl      (35,004 definition atoms, acyclic, 8,365 free vars)
    python3 -c "import harness" # -> orient.pkl  (bootstraps the orientation / topo order)
    python3 prop2.py           # -> prop2.json   (free=0 propagation, 38,998/39,033)
    python3 aeval.py ../best/new_instance_partial_39026.json   # atom view of the deliverable
Then: `engine.py` (cone eval + exact single-var solve), `fast.py` (incremental downstream
re-eval, 0.08 s/probe), `jclose2.py A B` (linear closure), `sparse.py`/`intsolve.py`/`js4.py`
(unit-pivot + HNF integer solve), `bitfeas.py B` (per-bit pin feasibility).

## Established results (all exact, reproducible)
1. Atom model is faithful: at the 39,026 deliverable exactly **8 atoms** are nonzero and they
   reproduce its failing set exactly.
2. **Seed all 8,365 free vars to 0 and propagate: only THREE atoms are violated** (38,998).
   With `x_18956 = C1` only three (39,009).
   C1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
3. The residual decodes exactly to: `OR(a,b)=1` **forced** (a=x_7715, b=x_34554, each an OR-tree
   over 178 / 78 free bits; x_9274=x_2300=1 is literal-pinned), plus a 2-way MUX (a20212) and
   a20215.
4. **The (1,1) branch closes all four core atoms analytically**: activate one free bit in each
   OR-tree, then set free vars `x_22162 = x_13682`, `x_30213 = x_18956 - x_32237` (4-iteration
   fixpoint).  With (x_4279, x_26005) the bad set becomes exactly the two bits' pin atoms
   [6668,12606,34497,34498].  (1,0) works too: `x_12186=x_13682`, `x_16742=x_18956-x_32237`.
5. Each activated bit carries pins `b*(free-K) = m*handle` and `pin = p*handle`
   (p = the 256-bit literal 115792089237316195423570985008687907853269984665640564039457584007908834671663 that appears throughout the file).  Single-variable greedy stalls at exactly 4 bad atoms.
6. **The repair is an exact linear Diophantine system.**  Closure around the pin set:
   4,008 vars x 2,996 atoms, linear except 1,810 (var,atom) pairs, built in 4 s.
   Small support: rationally feasible (rank(A)=rank([A|b])), **integer-infeasible**.
   The binding row for bit x_4279 is atom 19725: `-p*d_7815 = R` with `p | R` false, because
   atom 7128 forces `d_15951 = 0`.  So the obstruction is a p-divisibility (mod-p) condition,
   not a rank deficiency.

## Single highest-value next experiment
Run `bitfeas.py` over all 178 a-tree bits and 78 b-tree bits (it was mid-flight when stopped).
It answers, per bit, whether that bit's pin system is integrally solvable and names the binding
row.  Any bit whose system is feasible, paired with a feasible bit from the other tree plus the
(1,1) MUX assignment of item 4, gives a full witness.
