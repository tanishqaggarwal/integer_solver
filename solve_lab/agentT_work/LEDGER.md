# VERIFICATION LEDGER
**What this lab actually knows.**  Maintained by agent T (audit).  One row per load-bearing claim.

**How to read the `checked` column.**  This is the whole point of the document:

| mark | meaning |
|---|---|
| **T re-ran** | I executed it myself from cold and reproduced the number. Reproduction command given. |
| **T verified independently** | I established it by a *different route* than its author, usually F's certified-faithful parse. |
| **reported** | Recorded from the author's own file or from the coordinator. **I did not re-run it.** Treat as the author's claim, not as audited. |

Nothing below is stronger than its mark.  Where a claim exists in several atom numberings the row
says which — five are in play: **39,033** (F, K, T), **39,277** (P), **40,727** (E), **40,885** (I),
**42,267** (A, G, H), plus L's **9,032**-residual-atom engine.

---

## 0. THE RULES THAT EMERGED, each attached to the failure that produced it

1. **Check for repeats before reporting any rate.** — S, three times.
2. **The score counts equations, so price in equations.** — K, L, R, independently.
3. **A count derived from one parse is a fact about that parse until reconciled.** — Q, me, now M.
   My case: A's knob count is 9 at L=0 in a 42,267-atom parse and **24** in F's 39,033-atom parse,
   same window, same atoms, same variables.  The count moved 2.7x.  The 927 did not move under the
   same test, which is why it is a fact about the instance and the knob count is not.
4. **Validate a model by what it predicts is PRESENT, not by what it predicts is absent.** — R.
   Corollary I hit twice: an absence can be an artifact of your own convention (my 278) or of unit
   propagation's weakness (my 0/256 selectors).
5. **Dump the assignment and run the checker.** — mine, this round.  L's `|S|=2` closure was
   entirely model-internal until I dumped it; the step cost nothing and turned it into a fact.
6. **Never trust a symbolic expansion or a disjointness argument without direct recomputation.** — P.
7. *(mine, offered)* **Separate "this number is wrong" from "this result is wrong."**  O's exponent
   and atom count were both wrong and its Lemma was untouched; a wrong degree bound in L's solver
   cannot produce a false verified root.  Ask what the conclusion actually rests on before retracting.

---

## 1. VERIFIED

| claim | established by | checked | parse | reproduce | what would falsify it |
|---|---|---|---|---|---|
| **Deliverable = 39,026 / 39,033**, failing `[12231,12270,12350,14584,18673,22044,29125]` | prior campaign | **T re-ran** | raw | `python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json` | a checker-verified file scoring higher; none exists in the lab |
| **ker(M) = 0** — all 39,033 atoms forced | F | **T re-ran** | F 39,033 | `python3 agentF_work/peel_cert.py` | a nonzero atom vector in ker(M); the peel is a triangularity certificate, so only a parse error could |
| **Pivots are all ±1 or ±2** (37,889 ones, 1,144 twos) | asserted by F/FLEET | **T measured** — `peel_cert.py` only ever tested `pivot ≠ 0` | F 39,033 | `python3 agentT_work/t_pivots.py` | a pivot with an odd prime factor. Immaterial to ker(M)=0 over ℚ/ℤ |
| **M is faithful**: `eq_e = 0 ⟺ (Ma)_e = 0` | **nobody had tested this** | **T verified** — exact failing-set equality at 10 points incl. all-zeros, 4 partials, 2 random 30-digit | F 39,033 | `python3 agentT_work/t_faithful2.py` | one assignment where the atom-nonzero set and checker's failing set differ |
| **927 c>1 integer conditions is intrinsic** | L and P, unshared parses | **T verified independently, twice** — once borrowing only L's cofactor list, once borrowing nothing (13,092-atom family, 2.6x looser) | F 39,033 | `python3 agentT_work/t_927b.py` | the count moving under a further re-decomposition. It survived a 2.6x change; the knob count did not survive 2.7x |
| **The p-class**: exactly ONE atom in the instance carries the literal p, `(x26064 − p)`; its copy class is 220 wires and contains all six "shared slack factors" | Q raised it as open | **T proved it** — from `ker(M)=0` + faithfulness, every atom is zero in any solution, so the pin propagates | F 39,033 | `python3 agentT_work/t_slack3.py` | a second literal-p atom, or a break in the copy chain |
| **The coordinate hand-off is a one-atom affine alias on ALL 764 parent/child links, slack = p·u in 764/764** | Q (573 links) | **T verified over all 764** after correcting my own pairing | F 39,033 + L's calib2 | `python3 agentT_work/t_cross.py` | a link with no one-atom alias, or slack that is not p·u |
| **Cancellation is a VALUE property, not a support property** | L | **T verified from the opposite direction** — from the deliverable, not through L's constructor: zeroing 12 cofactors keeps the atom support byte-identical and moves 7 → 12 | F 39,033 | `python3 agentT_work/t_cancel.py` | identical support ever forcing identical cost |
| **Cofactor criterion** `e contains a ⟺ u_a ∈ vars(e)` | L | **T verified exhaustively**: 3,681/3,681 free, 3,681/3,681 in exactly one atom, 3,681/3,681 with `eqs(u) == eqs(atom_u)`. Zero violations | F 39,033 vs checker varsets | `python3 agentT_work/t_cofactor.py` | any cofactor occurring in two atoms |
| **L's `\|S\|=2` closure over ℤ** — all 927 discharged; **39,018/39,033**, exactly 2 nonzero atoms (the target congruences), their 15-equation footprint **equal to** checker's failing set | L (model-internal only) | **T reproduced, dumped and checker-verified** — L's script dumps nothing; `assign_L2.json` predates the run by 2h | L 9,032 → F 39,033 | `python3 agentT_work/t_S2.py && python3 agentT_work/t_S2b.py` | a nonzero atom outside the two congruences, or a failing equation outside their footprint |
| **O's Lemma: `S⁴ = 0 ⟹ S = 0`, unconditional** — S affine in all 43 vars, `dS/dx_4432=+1`, `dS/dx_28730=−1`, a23618 at coefficient +1 | O (as "T²", 20 atoms) | **T verified, with 2 corrections**: the equation is **S⁴ not S²**, the form has **18 atoms not 20**. F's parse independently gives the same 18 coefficients | raw + F 39,033 | `python3 agentT_work/t_eq8680.py` | a factor that is not a repeat, or S non-affine in some variable |
| **N's detach exhaustion: the 2^65 detach lattice has exactly 16 states, OPT = 5 for all 16, best 39,026** — because only 4 of 65 pool vars have witness ≠ gate (`{642,28730,29854,31864}`) | N | **T verified independently** — in F's parse a pool var's defining atom is nonzero iff witness ≠ gate; exactly those 4 are nonzero. **And T closed the gap N left**: 0 of the other 61 reach a witness var anywhere in the full 30,001-definition DAG, so the no-op holds at all 16 states, not just the witness. **16 by proof; signatures complete by construction** | F 39,033 | `python3 agentT_work/t_detach.py` | a 62nd pool var with a nonzero defining atom, or one of the 61 depending on a witness var |
| **M's exhaustive placement enumeration**: all 4,096 subsets of `H12` priced, nothing above **39,026**, witness unique at support 4; all 30 subsets attaining 39,026 are supersets of the witness (M's `verifysup.py`) | M | **T verified the ENGINE independently** — driven on the witness subset it reproduces 39,026 with the deliverable's exact 7 failures under `checker.py`, F's parse gives the deliverable's 7 atoms, and the assignment is **byte-identical to the deliverable (0 of 38,748 vars differ)**. Exact at **9/9** spot-checked subsets spanning 39,008–39,026. **Scope: 9 of 4,096 checked — the verdict on the other 4,087 rests on the scorer, not on my check** | M's engine → checker + F 39,033 | `python3 agentT_work/t_meng.py && python3 agentT_work/t_meng2.py` | M's engine disagreeing with `checker.py` on any subset |
| **≥7 failing is unconditional and exhaustive at L=0** over all 24 genuine knobs (A used 9); weight ≤6 admits nothing, weight 7's unique lightest set **is** the deliverable's | **T** (new result) | self-tested end-to-end by that coincidence | F 39,033 | `python3 agentT_work/t_conda.py 0 8` | a weight-≤6 mod-p-consistent violated set |
| **Q's ladder**: 249/249 checkable doublings exact, 253 distinct points, the 3 "missing" exponents are decoded leaves (x18184, x22579, x33434 = 2⁵¹G, 2¹⁷⁶G, 2⁴¹G), `N·G = O` | Q | **T re-ran** — stronger than Q claimed; nothing is inferred | Q's model | `python3 agentT_work/t_ladder.py` | a failed doubling, a duplicate point, or an exponent with no decoded leaf |
| **All 256 leaf selectors are booleans and free variables** | Q/L | **T re-ran** | E's parse | `t_fold2.py` | a selector with a definition |
| **No slot output feeds another slot directly: 0 of 764** | Q (0/383) | **T verified independently** from L's calibrated model | L's model | `python3 agentT_work/t_alias.py` | one direct link |
| **H's 722 dormant handles are a safe exclusion** — census identical at 8 configurations incl. 3 random 8-selector ones | H (measured at 1 state) | **T re-measured at 8** | H/E's parse | `python3 agentT_work/t_dormant2.py` | one handle changing class at some configuration |
| **I's eq8680 hunt conclusion** — the 5 nonzero-effect groups I never tested all give minfail > 6 | I (coverage claim was wrong, see §3) | **T re-ran using I's own solver** | I 40,885 | `python3 agentT_work/t_hunt_gap.py` | a group with minfail ≤ 6 |
| **A's linearity filter discards nothing** — `knobs_raw == knobs_linear` at every level 0–6 | A | **T re-ran A's own `eqwin.py`** | A 42,267 | see RESUME_T §A6 | a level where the two differ |

---

## 2. CONDITIONAL — scope is in the row, not a footnote

| claim | scope it is true in | established / checked | what would move it |
|---|---|---|---|
| **M's useful enumeration space is 2¹⁸** | after **T's correction from 15 → 18**: L's "exactly 15 incident of 3,681" was applied to a census that omits 33 genuine p-handles whose guards are stage checks and leaf pins, **3 of which are incident** (x10422, x15120, x35531 — all three appear as terms of O's `S`) | L; **T corrected** (`t_final.py`) | re-running the incidence filter over the full p-handle family (3,707, or 3,714 both operand orders). **Note: `S` also has 18 terms. Different 18s — do not conflate** |
| **L's cancellation freedom is 4-dimensional** (`x1329 +3, x9413 +4, x10903 +3, x17325 +4`) | 8 of L's 12 named cofactors are **already 0 in the deliverable**, so they are no-ops. The gap is **5 (7→12)**, not 6, and the far side is **12**, not L's 13 | L; **T corrected** | note the 4 h-wires `x642, x28730, x31864, x29854` are *also* effectively assignable because the deliverable already breaks their defining atoms — **true dimension is > 4 and unestablished** |
| **The reduction closes** | **MOD P only.** The hand-off slack is `p·u`: identically zero mod p, nonzero over ℤ. Over ℤ what remains is exactly the 927 | Q, L; **T supplied the mechanism** | any statement of the reduction claiming closure **over ℤ** without discharging the 927 is wrong as written. L's says mod p and is fine |
| **Q's six restored sweeps** | valid because the fold is a group operation **mod p**, and the p-slack vanishes there — so the group model is exact in the ring the sweeps work in | Q; **reported**, T did not re-run the sweeps | T *did* verify the ladder and the group order they rest on (§1) |
| **O's seven-way 1-for-1 trade** | over its **34 inputs**, **frame B's orientation** | O; **reported — T's audit is PARTIAL and reached no verdict.** T found (a) `K` rebuilt in the default orientation gives **23**, not 34 — frame-dependent, flagged not refuted; (b) the uniformity is **not** structurally forced (7 knobs move a failing row with `dS=0`), so it is a genuine search result resting on O's collateral accounting, which is **unaudited** | reproducing frame B and checking the collateral accounting |
| **K §4** | **premise unrefuted, NOT established.** K's own file carries a DO NOT REBUILD THE GUARD header | K; **reported, T did not re-run** | — |
| **N's OPT = 5 / `outside = 0` pricing, and the 924/924 p-obstruction** | under **`fwd2`'s orientation** and the witness region (|R| = 12/13). **T did not re-run these** — only the 16-state reduction they sit on. Independent of O's Lemma: the region excludes eq8680, and N established that O's Lemma is instead exactly the 39,025→39,026 step | N; **reported** | N's own re-orientation run |
| **L's degree ≤ 3 bound** | real, not a 5-point aliasing artifact — **T re-fitted at 7, 9, 11 points, same top degree every wire**. But it bounds **cost, not correctness**: the recomputation guard rejects a bad root, so a wrong bound can only cause a *missed* solution | P, L; **T verified** (`t_deg.py`) | nothing about it can invalidate a verified result |
| **L's closure generalises** | **ANSWERED — NO, and the cause is localised.** T ran the sweep on a nested chain: `|S|` = 2 (control), 3, 5, **6, 7** all close (2 nonzero atoms, **39,018**, identical 15-equation failing set); **`|S|`=8 fails** (3 atoms, **39,002**). **Not a size horizon — a single leaf, `x34974`.** Its residue `((x21408*x10138)-(15333171*x658))`, c = 3·7·19·83·463, has **no root on any of its 6 candidate wires and is blocked by collateral on none** — so it is **NOT** the `\|S\|`=17 shared-wire simultaneity but is consistent with **L's bivariate residue**. Scope: one ON-set per size; single-wire granularity | T ran it; `python3 agentT_work/t_sweep2.py && python3 agentT_work/t_leaf.py` | a two-wire shift clearing the residue, or another 8-leaf set that closes |

---

## 3. WITHDRAWN — do not re-derive

| claim | withdrawn by | why |
|---|---|---|
| **The five barriers** | their authors | each was a property of a filtered knob set or an unstated configuration reported as a property of the instance |
| **Q's six sweeps** | Q, then **restored** | withdrawn on a term-count mismatch; restored after Q reconciled 47,198 terms against F's 39,033. Now conditional (§2) |
| **K's partition-theorem withdrawal** | K, then **retracted the withdrawal** | K's null searched for a *direct* composition; **T showed 0 of 764 links are direct**, so the null was consistent with the alias layer. K then found the alias inert in that assignment and the match literal — so the withdrawal was wrong for a *different* reason than T proposed. Net: theorem stands, T's alias explanation superseded |
| **S §3** ("a20215 ≢ 0 is unreachable") and **S §6i** | S | S's own scope limit. T's contribution: refuted the "different spaces" hypothesis — Q's 256 leaf selectors and S's 256 cluster booleans are the **identical set** |
| **R's accumulator model** | R | **reported**; the rule it produced is rule 4 above |
| **M's 32-handle pool** | M | **reported** — superseded by M's own H12/H16 incident sets |
| **L's "15 incident atoms"** | superseded by **T** | census omitted stage-check/leaf-pin guarded p-handles; true count 18 (§2) |
| **T's "278 multi-hop aliases"** | **T** | my own artifact. I matched `OUT[n][j]` to `OUT[child][j]`; L's `calib2` had already measured the coordinate alignment (188 orient=1, 67 orient=0). With the cross allowed: 486 + 278 = 764, **0 unaliased** |
| **T's "liveness is not determined by the selectors"** | **T** | Q's slot analysis: the pin is `sel*(w−C) − z`, so routing *is* determined — by a simultaneous system, not propagation. The narrower **"forward evaluation from the selectors does not realise an ON-set"** stands (0/256 at four bases) |
| **A's THEOREM.md §7** ("any assignment satisfying all equations makes every atom zero") **as applied to A's own atoms** | **T** | true in F's 39,033 parse; **false** in A's 42,267 (kernel dim ≥ 3,234) and I's 40,885 (≥ 1,852). Rule 3 |

---

## 4. GENUINELY OPEN

1. **L's bivariate residue** — the last undischarged condition at `\|S\|=17`, and now the likeliest reading of T's `\|S\|=8` residue too (no univariate root on any of 6 wires). **Concrete next test: 15 wire-pairs, bivariate root-find.** The system is nonlinear (a `p·t_w·t_v` term survives mod c), so no linear solve expresses it.
2. **The collision criterion.**
3. **The scalar recovery itself** — a 256-bit discrete log in a prime-order group with no exploitable structure. Q's position, which T's ladder check supports: ~2^128, and no amount of circuit decoding reduces it.
4. **What the 278 crossed-index aliases' partner wires route through** was resolved (§3); but *whether every one of the 764 aliases is inert over ℤ* is exactly the 927 question (§2).
5. **Whether tuning the 4 real cofactors beats 7** — M's search. True dimension unestablished (§2).

---

## 5. POINTERS
Full reasoning for every T row is in `agentT_work/RESUME_T.md` (740 lines, sections A–AE); it is
organised by pass, and each section names its scripts.  All T scripts are in `agentT_work/` and
were run read-only against other agents' directories with `PYTHONDONTWRITEBYTECODE=1`.
