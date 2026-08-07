# RESUME_E — agent E private checkpoint

**WITHDRAWN ANGLE:** generator/authorship reverse-engineering.  Do NOT resume it.
Everything below is computed from the equations as mathematical objects.

## Verified
- `best/new_instance_partial_39026.json` = 39026/39033.  CONFIRMED (18 s).
- My own best so far: **39,015** (fails=18).  Nothing of mine reaches 39,026 yet.

## Machinery (agentE_work/)
- `parse3.py`->`model3.pkl`: 39,033 eqs = outer * Z-combo of 40,727 atoms (9,710 eqs are S*S).
- `dag.py`->`dag.pkl`: 35,004 atoms are defs `x_out - RHS`, acyclic; 8,365 free vars.
- `harness.py`/`engine.py`: forward free-inputs -> all vars (0.15 s); cone eval (0.01 s);
  `E.solve_for(atom,var,seed)` exact single-var solve.  Orientation cached in `orient.pkl`.
- `fast.py`: **incremental** downstream-only re-evaluation, 0.08 s/probe, verified exact.
- `jclose2.py A B`: closes the exact linear map free-vars -> atom residuals (4 s).
- `sparse.py` + `js4.py`: unit-pivot elimination then HNF (`intsolve.py`) integer solve.

## Established
1. **free=0 => only 3 violated atoms** (38,998).  `x_18956=C1` -> 3 atoms, 39,009.
   C1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
2. Core residual = `OR(a,b)=1` forced (a=x_7715, b=x_34554; OR-trees over 178 / 78 free bits)
   + a 2-way MUX (a20212) + a20215.
3. **(1,1) branch closes all four core atoms exactly**: set one free bit in each OR-tree, then
   `x_22162 = x_13682`, `x_30213 = x_18956 - x_32237` (4-iter fixpoint).  With bits
   (x_4279, x_26005) the bad set is exactly the two bits' pin atoms [6668,12606,34497,34498].
   (1,0) branch: `x_12186 = x_13682`, `x_16742 = x_18956 - x_32237`; costs 4 atoms too.
4. Pin gadgets `b*(free-K) = m*handle`, `pin = p*handle` (p = secp256k1 prime).  Single-move
   greedy stalls at exactly 4 bad atoms; residual relocates down chains (`runs/chain_A.log`).
5. **Closure around bad=[6668,12606,34497,34498]: 4,008 vars x 2,996 atoms, exactly linear
   except 1,810 (var,atom) pairs.**  round<=1 (61 vars): rationally feasible, integer-INfeasible.
   round<=2 (342 vars, 335 rows): rationally feasible (rank 297=297); integer solve running.
   => the small-support obstruction is INTEGRALITY, not rank.  Grow support to defeat it.

## Next experiment (single highest value)
`nohup python3 js4.py jac_4279_26005.pkl 1 5 js4_A.json` — sparse integer solve at growing
support.  If a round is integer-feasible, apply and re-verify with `../checker.py`.
Then sweep bit pairs: `jclose2.py A B` + `js4.py jac_A_B.pkl 1 5 out.json`;
per-bit costs in `runs/search2.log`, a-bits and b-bits listed there.
