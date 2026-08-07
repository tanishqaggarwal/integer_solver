# Agent H log — integer/polynomial analysis of EQUATIONS.txt

## Step 0
Read PROMPT.txt, RESUME.md, STATE.json.  Verified the claimed best:
`python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
-> satisfied 39026/39033, failing [12231,12270,12350,14584,18673,22044,29125].  CONFIRMED.

## Step 1 — the frame (kept; this is the durable asset)
Independent parse (`model.py`): 42,267 atoms / 39,033 equations.  Orienting every atom of
syntactic form `x_t - rest` as the definition `x_t := rest` gives an ACYCLIC gate DAG over all
38,748 variables: 8,747 free inputs, 30,001 defined, 12,266 check atoms, ZERO cycles.
(The lab's prior frame had 7,273 free inputs and 1,800 variables inside gate cycles.)
Forward evaluation from ALL-ZERO free inputs scores **39,005** with only **5** nonzero check
atoms — checker-verified on `allzero_fwd.json`.

## Step 2 — decomposition hypothesis, REFUTED by my own measurement
eq-var bipartite graph: 1 component.  atom-eq graph: 1 giant + 3,234 singletons.
Free-input hypergraph (hyperedge = equation): 1 component, 8,747 vars / 18,248 equations;
20,785 equations are identically satisfied by forward evaluation.
Residual closure = 6,007 free inputs / 6,026 checks / 9,244 equations.
=> NO separator, NO block decomposition.  My assigned angle's central hypothesis is FALSE.

## Step 3 — the residual, in raw terms
Exactly three conditions remain:
  (i)  x_9274 = 1   (x_9274 = (x_7715 + x_34554) - x_7715*x_34554, over 256 selector variables)
  (ii) x_37892 = C1 (mod p)
  (iii)x_13682 = C2 (mod p)
with p = 115792089237316195423570985008687907853269984665640564039457584007908834671663,
C1 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626,
C2 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002.
Forced because a688 = 8863713*(x_18956 - C1) - x_14257 with x_14257 = p*x_7497 and gcd(8863713,p)=1;
a1618 = x_24468 - C2 - p*x_11436; a30982 = x_18956 - x_37892 - p*x_22820;
a30980 = x_24468 - x_13682 - 12354891*p*x_14393.

## Step 4 — how the 512 load-pin constants enter (measured, not assumed)
Each selector variable b carries exactly two pin atoms `b*(x_F - H) - c*x_handle` (512 constants,
`leafpins.json`).  They are combined by a tree of mux atoms.  For two active selectors the
combined pair (X3,Y3) is forced by the polynomial identities
   (X1+X2+X3+K)*(X2-X1)^2 = (Y2-Y1)^2 ,  (Y3+Y1)*(X2-X1) = (Y2-Y1)*(X1-X3)   (mod p)
K = 97553848499418123410591666447050222001188385549510401465815187079080512838891.
Measured 56/56 on (u,w) selector pairs, with a 1-off control failing 56/56.
Over Z the triple additionally needs 2264251 | (x_15286/p), 2264251 = 11*43*4787; solving that
congruence for one lift parameter closes all three checks EXACTLY over Z (3/3).
**Consequence: the constants combine NONLINEARLY.  No linear relation among them expresses the
residual, so linear/lattice relation-finding over the constant table cannot be the operative tool.**

## Step 5 — tools built (kept)
`fast.py` incremental exact evaluator, 1.4 ms per move.
`close2.py` bottom-up cascade closer: closes the entire pin tree in ~32 exact assignments with NO
search (smallest-atom-support first, freeze each assigned variable).
`close3.py` adds an indirect phase.  `scan1.py`: every single-selector state closes to **39,018**
with exactly 4 nonzero atoms {30980, 30982, 36185, 40812}.
Beam search over repairs (`beam2.py`) plateaus at 39,012 and oscillates — the constructive closer
strictly dominates it.
