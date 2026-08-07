# RESUME_E — agent E checkpoint

**FRAMING (per coordinator):** no curve / point / group / ladder / discrete-log vocabulary is
used anywhere below or in LOG.md.  `p` denotes only the 256-bit integer literal
115792089237316195423570985008687907853269984665640564039457584007908834671663 that occurs
verbatim in EQUATIONS.txt.  Every claim is about integer congruences and integer-linear
relations among the given polynomials.  (Also previously withdrawn and NOT resumed:
generator/authorship reverse-engineering.)

## Verified scores
- `../best/new_instance_partial_39026.json` = **39026/39033** (fails
  [12231,12270,12350,14584,18673,22044,29125]).  CONFIRMED with `../checker.py`.
- Mine: 106 checker-verifiable states at **39017/39033** (`bitsol_<bit>_39017.json`, spot-checked
  with `../checker.py`).  Best of mine overall: 39,015-39,017.  Baseline not beaten.

## Rebuild (about a minute)
    python3 parse3.py; python3 dag.py; python3 -c "import harness"; python3 prop2.py
Modules: `engine.py` (forward map + cone eval + exact single-var solve), `fast.py`
(incremental downstream-only re-evaluation, 0.08 s/probe), `sparse.py`+`intsolve.py`
(unit-pivot elimination + HNF integer solve), `bitfeas2.py` (per-bit closure+solve),
`scanfork.py` (tree sweep, hard per-bit kill), `iterfix.py` (iterated closure+solve),
`subsetfeas.py`, `big01.py`.

## Structure (exact, reproducible)
1. 39,033 eqs = outer scalar x Z-combination of 40,727 **atoms**; 9,710 eqs are squares S*S.
   35,004 atoms are definitions `x_out - RHS`; the definition graph is ACYCLIC; 8,365 free vars.
2. **Seed all free vars to 0 and propagate: only THREE atoms violated** (38,998/39,033).
3. The residual decodes to: `OR(a,b)=1` forced (`a=x_7715`, `b=x_34554`, each an OR-tree over
   178 / 78 free bits), a 2-way MUX (a20212), and a20215.
4. Activating one free bit costs exactly 2-3 "pin" atoms of the form
   `bit*(free - K) = m*handle` and `pin = p*handle`.
5. **Pin repair is an exact linear Diophantine system**; the obstruction is p-divisibility
   (integrality), not rank.

## THE SCAN (complete)
Per bit: closure of free-vars -> atom residuals, then exact integer solve; every solution
re-verified by exact evaluation.
* **a-tree, 178/178 bits: 56 FEASIBLE**, 102 core-infeasible, 6 infeasible on an explicit
  p-divisibility row (e.g. `x_1746 -> row 33794: rhs % -p != 0`), 8 undecided (HNF core >
  guard), 6 timed out.  (`runs/scanA6.log`, `scanfork_A.pkl`)
* **b-tree, 78/78 bits: 50 FEASIBLE**, 28 infeasible (26 re-confirmed with the improved
  solver).  (`runs/scanB.log`, `runs/recheckB.log`, `scan_B.pkl`)
* **All 106 feasible bits verify exactly with ZERO residual atoms outside the selector core**
  (`exact=(16, [])` -> 39017/39033).

## SUBSETS — the answer (SUPERSEDED IN PART, see "THE TRIPLE IS SOLVABLE" below)
| bits on | residual atoms after iterated exact repair |
|---|---|
| 1 | **0, 0, 0** |
| 2 | 3, 3, 4 | 3 | 7, 7, 7 | 4 | 6, 7, 9 | 5 | 9, 10, 10 |

Every subset of size >= 2 turns on a product ("AND") node at the meet of the chosen bits in the
OR-tree and leaves an irreducible triple, e.g. for (x_1530, x_1603) with `x_24195 = 1`:

    x_24195*(5002401*U + 15322661*V)           = 0            (exact over Z)
    x_24195*(15944455*U + 4826103*V)           = p*x_34496    (congruence mod p)
    7952523*x_24195*(14913407*U + 11707765*V)  = p*x_3193     (congruence mod p)

with `U = x_29210 = x_25848 - x_17317`, `V = x_8736 = x_18682 - x_28841`.  The 2x2 determinant
of the two mod-p forms is 15944455*11707765 - 4826103*14913407 = 114700293930154, coprime to p,
so the pair has rank 2 and forces `U = V = 0 (mod p)`.  Iterated closure+solve is a fixed point
on these triples.

**Independent or only through the sum?  NEITHER.**
* Two individually-feasible bits move 3 variables in common — `x_14853, x_31339, x_6083` — and
  require **conflicting** values on all three (0 of 3 agree).  Across the 50 feasible b-bits the
  50 required values of `x_14853` are pairwise distinct, and distinct mod p.
* The coupling is not additive: with both bits on, `U` and `V` are not the sum of their
  singleton values, over Z or mod p; the product flag `x_24195` flips 0 -> 1.
So single-bit solutions cannot be recombined linearly, and each subset poses a fresh rank-2
congruence system.  This is a measurement over the supports searched, not a proof.

## Highest-value next experiments
1. The 8 "core too large" + 6 timed-out a-bits are the only undecided singletons — raise the HNF
   budget (`sparse.solve_sparse(..., maxcore=)`) or use a sparser integer solver on them.
2. Attack the triple directly: solve `U = V = 0 (mod p)` together with
   `5002401*U + 15322661*V = 0` over Z, using the free variables in the cones of
   `x_25848, x_17317, x_18682, x_28841` rather than the selector bits.  `runs/iterpair.log`
   shows the iterated *linear* model stalls there; a targeted 2-unknown congruence solve has
   not been tried.

## THE TRIPLE IS SOLVABLE — subset barrier REFUTED (supersedes the subset verdict above)
The size>=2 obstruction triple closes exactly.  `p|U` + the exact row are the only two
independent conditions (p|V follows), and `(U,V)` is exactly affine in the two non-boolean
knobs **x_30468, x_33169** — which every earlier search missed because they lived in rows the
closure dropped as nonlinear.  `triple4.py` solves them in closed form (one free parameter k,
second congruence has modulus 1 so every k works); the two quotient handles `x_34496`,
`x_3193` then kill atoms 722/724.  `triple8.py` then solves **44 of 46** rows of the affine
system simultaneously -> `triple8_39005.json`, **39,005/39,033 with only TWO nonzero atoms
in the whole instance** (verify: `python3 verifyE.py triple8_39005.json`; the values exceed
Python's 4,300-digit string cap, which is all verifyE raises — it calls checker.py unmodified).

**The instance now binds in exactly two rows:**
    a20215 : x_24530 - x_5647*x_24908           ->  x_24908 = C1
    a28647 : x_36433 - (x_36990 + x_19239)      ->  x_26386*x_6083 + x_27475*x_33708 = C2
C1 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
C2 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
Both say "an accumulator must take one specific 296-bit value".  They also blocked the
single-bit branch (LOG.md §11), so they are frame-independent within my model.
NEXT: close a20215 / a28647 by an affine solve over the free variables in the cones of
x_24908 and of (x_6083, x_33708) — never attempted; the knob search there must include the
non-boolean integer handles, which is exactly what unlocked the triple.

## a20215 / a28647 (the two rows the instance now binds on) — RUN, and the reason it stops
The cluster closes into an exact 5-row system over 8 integer knobs (LOG.md §14).  Infeasible,
and the mechanism is measured: `c1 = d(a10187)/dx_37012` and `c2 = d(a20212)/dx_14393` are both
**= 0 (mod p)**, so a10187 pins `d_31339 = 0 (mod p)` and a20212 pins `d_14853 = 0 (mod p)`,
while a20215 and a28647 need those knobs to hit `-R1` and `R2` mod p, and
R1 mod p = 22981624690591324143788809642515852940280603493270692712106986169263210356252,
R2 mod p = 44159679639019146557987083382852396884224992023970032213706899677695745279353,
both non-zero.  Every knob reaching a10187/a20212 enters mod-p trivially; the only non-boolean
movers of a20212 are x_11436 (coefficient exactly p), x_14393 (c2 = 0 mod p) and x_14853 itself.

**NEXT:** find a knob reaching a10187 or a20212 with coefficient COPRIME to p.  That is the
single quantity that decides this cluster.  Search beyond the cones tried (326 candidates so
far), and treat the 181 boolean selector movers of a20212 as 0/1 decisions rather than affine
variables — the boolean class was excluded by the affineness filter, and excluding a knob class
by filter is exactly the error that produced the earlier false barrier.

Verification rule for my states: values exceed Python's 4,300-digit string cap, so
`checker.py` cannot PARSE them.  Use `python3 verifyE.py <file>` — it raises only that cap and
calls checker.load_equations / load_assignment / evaluate_all unmodified.  Say so explicitly
in any report; a bare "checker.py says" would be false for these values.

## CORRECTION to the a10187/a20212 verdict — the boolean class was excluded, and it matters
Measured by exact re-propagation at 0/1 (`boolknob.py`): **23 boolean flips have a residual
delta nonzero mod p on a10187 and 178 on a20212.**  My §14 statement that every knob there
enters with a coefficient divisible by p is FALSE — it held for the affine class only.
The cluster is NOT mod-p sealed.
What still blocks: the nonzero deltas carry only 1 distinct residue on a20212 and 3 on a10187,
while the required multiplicities are full 256-bit numbers (k, m in LOG §15.1), so counting
flips cannot supply the shift; and each flip drags in its own pin atoms, which the 8-knob
cluster solve cannot reach (`boolsolve.py`), while flip + full pin repair tops out at 39,005
(`boolfix.py`).
**NEXT:** compose properly — take a bit from the 106 whose pin system is known integrally
solvable (`bitsol_*.json`, `scan_B.pkl`, `scanfork_A.pkl`), apply its pin solution AND the
cluster's 8 affine knobs in ONE simultaneous system rather than sequentially; the sequential
composition is what fails, not either half.  Also: probe whether flipping a bit changes the
three a10187 residues (they are the scarce resource), since §15.2 shows the residues are
context-dependent.
