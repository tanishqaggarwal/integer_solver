# RESUME_S — agent S (lattice methods), self-contained handoff

## 0. Verification rules
- Baseline `../best/new_instance_partial_39026.json` = **39,026/39,033** re-verified by me with
  plain `solve_lab/checker.py` (failing `[12231,12270,12350,14584,18673,22044,29125]`).
- States with >4,300-digit integers must be verified with `../agentE_work/verifyE.py`, NOT
  `checker.py` (which raises ValueError on parse). Say which was used.
- **Every "cannot move X" claim below states its knob set AND selector configuration.**

## 1. Setup (rebuild in ~5 s)
Work from `solve_lab/agentS_work/`; symlinks to E's `orient.pkl`, `users.pkl`, `triple8_seed.json`
are already here. `common.py` loads agent E's forward engine read-only.
    python3 -c "import sys;sys.path.insert(0,'.');import common"   # 3 s
Scripts: `meas.py` (delta table, unfiltered), `table.py` (exact supports, -> table_cfg0.pkl),
`hcheck.py` (per-atom pure-handle search), `hrepair.py` (handle-repaired scoring),
`local1-3.py` (deliverable's local system), `cascade.py` (exact cascade repair, eq-level score),
`sat.py` (saturation test), `image.py` (selector->residual image, -> image.pkl),
`lattice.py` / `lat2.py` / `lat3.py` (the exact integer lattice solve).

## 2. THE MAIN RESULT — the endgame condition, computed exactly
At cfg0 (= `triple8_seed.json`, selectors x_1530=x_1603=1; residual atoms a20215, a28647):

Build the exact integer system `D n = -R` over **54 affine knobs** — the complete affine knob set
of the 5-row cluster cone, *including every knob that moves only one row* (x_22820, x_11436,
x_14393, x_26489, x_37012, x_18956, x_31339, x_14853, x_6083, x_30468, x_33169), plus the pure
handle of every atom in play, plus the cone free vars of the no-handle atoms. 49 atoms in play,
only **two have no handle: a726 and a28647**.

`lat3.py` result: eliminate a20215/a28647, solve the other 47 rows exactly over Z (FEASIBLE,
kernel dim 7), then compute the image lattice on (a20215, a28647).
**The reachable lattice is exactly p·Z^2**, p = 2^256 - 2^32 - 977. Membership of the required
offset: **NO**. So at cfg0 the endgame condition is exactly

        a20215 == 0 (mod p)   AND   a28647 == 0 (mod p)

and both are nonzero: 22981624690591324143788809642515852940280603493270692712106986169263210356252
and 44159679639019146557987083382852396884224992023970032213706899677695745279353.
This is a *derived* lattice statement over the unfiltered knob set, not a filtered barrier.

## 3. Why the handles are the whole story (new, exact)
1,815 free variables are **pure handles**: they move exactly one atom, affinely, by exactly ±p
(`handles2.py`). Every atom in play except a726 and a28647 has one. So an atom is satisfiable iff
its residual is ≡0 modulo its handle step. Measured steps at cfg0:
a20215 p, a20212 p (and 12354891p), a7389 p, a10187 **7942211*p**, a747 p, a32148 p, a28035 p,
a28037 p, a20652 6672769*p, a28033 8640431*p.  a28647, a30787, a26958, a40306, a726: **no handle**.
The trade knobs are exactly 1-for-1: x_14853 (a20212 -1, a28647 +1), x_6083 (a7389 +1, a28647 -1),
x_31339 (a10187 +1, a20215 -1), x_18956 (a747 +8863713, a20215 +1).

## 4. SATURATION — reproduced independently, and it is STRONGER than reported
`sat.py`, measuring the atom values directly (not bad-atom-dict deltas): at cfg0 the 219 moving
booleans fall into 2 classes (178 / 41) on (a20215 mod p, a28647 mod p), and
**every tested subset — within-class AND CROSS-class — gives the delta of a single member.**
Cross-class saturation was not previously recorded. Consequence: the selector-to-residual map is
a priority/one-hot selection, not a sum, so it is NOT a subset-sum and LLL on the bit vector is
the wrong tool (this confirms, by a second route, why the earlier LLL attempt failed).

## 5. The image of the selector map is TINY (`image.py`, image.pkl)
400 random configurations over all 256 cluster booleans (random on-sets of size 1..200, existing
ON bits randomly cleared): only **18 distinct 5-tuples** of (a7389,a10187,a20212,a20215,a28647)
mod p. a20215 mod p takes 2 values, a28647 mod p takes 2 values, **never 0**, and the combined
invariant (a28647 + a20212 + a7389) mod p takes 11 values, never 0.

## 6. The deliverable is a CANCELLATION CODE — its 7 is not 7 broken atoms
`local1.py`/`fix1.py`: the 39,026 assignment has 8 nonzero atoms whose footprints union to **12**
equations; 5 of those 12 are satisfied by cancellation. Removing any subset breaks the code:
setting x_31864=0, x_10903=0 (which zeroes a36661 and a36663 with provably zero side effects,
occ=2 and occ=1) makes it **worse: 11 fails / 39,022**. Depth-4 exact cascade repair from the
deliverable (`cascade.py`) found nothing better than 39,026.
Objective note: **769 atoms have equation-footprint exactly 1, covering 769 distinct equations** —
so a state whose only nonzero atoms are <=6 footprint-1 atoms would score >=39,027.
The deliverable's residual atoms are all `x_A - p*x_private` handle atoms with x_17499 = x_22665 =
x_28599 = x_28961 = p exactly; a36663 is `x_31864` and a36662 is `x_7075 * x_8731`.

## 7. Scores
- Best verified by me: **39,026** (the existing deliverable, re-verified with `checker.py`).
- No improvement found yet. E's 39,005 reproduced exactly.

## 8. Next experiments, priority order
1. **`lat4.py` OOMed (exit 137)** running flint HNF on 5,000-bit entries across 18 image points.
   Re-run with the residuals REDUCED MOD p first (the lattice is p*Z^2, so all that matters is the
   mod-p class) — that makes every entry <=256 bits and the HNF trivial. Highest EV.
2. **The two non-affine knobs x_30163 and x_11559 are the only unexplored freedom.**
   x_30163 moves a28647 (296 bits, NOT a multiple of p) and a30787; x_11559 moves a10187, a26958,
   a40306. a30787, a26958, a40306 have NO handle. Sweep their values and see whether a28647 mod p
   moves. If it does, section 2's barrier is broken. THIS IS THE ONE TO DO FIRST after (1).
3. Re-run the section-2 lattice analysis at each of the 18 image points (configuration-dependent:
   both the handles and the targets move — re-measure, never carry over).
4. a726 is the other no-handle atom; check whether it is genuinely forced or just unhandled.
