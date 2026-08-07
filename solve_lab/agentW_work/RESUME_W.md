# RESUME_W — agent W.  Round 2: THE CLASSIFICATION QUESTION.  Answered.

**Best verified score: 39,026 / 39,033** — `solve_lab/best/new_instance_partial_39026.json`,
re-verified from cold this session with `solve_lab/checker.py`, failing
`[12231, 12270, 12350, 14584, 18673, 22044, 29125]`.  **I did not beat it.**

Environment rebuilt after the restart: `python3 model.py; python3 fwd2.py` in `agentW_work`.
`model.pkl` md5 `ca31583ff6d604bfdc5b72da2f0e3f84`, `fwd2.pkl` md5 `69a85eae6b52612be2bbabc2aae2f8f6`
— **md5-identical to agent N's**, so this is H's model, not a re-derivation.  `agentH_work`
untouched (it still has no `.pkl`).  `PYTHONDONTWRITEBYTECODE=1`.  No git commands.

Round-1 results (frame-B budget, |K|=34, the 32-way trade, integrality) are unchanged and are
in check-in 92; the round-1 section of this file is preserved at the bottom.

---

# THE ANSWER: the classification closes at TWO families.

## 0.  The ring, stated first

I did **not** assume the mod-P picture.  I measured the integer form of every block.
Every one of the 383 law blocks carries **five** atoms, not three:

```
congruence k=1,2,3 :   a_k  * L     * ( c_k1*N1 + c_k2*N2 )  =  c_k  * P * u_k
off-pin    j=5,6   :   a'_j * (1-L) *   i_j                  =  c'_j * P * u'_j
```

with `N1 = E*A^2 - B^2`, `N2 = A*(i3+i6) - B*(i2-i5)`, `A=i1-i2`, `B=i4-i3`, `E=i1+i2+i5+Q`,
`L` the block's liveness gate, and `u_k, u'_j` **private handle variables**.
`P = 115792089237316195423570985008687907853269984665640564039457584007908834671663`
(= `x_26064`, 220 aliases), `Q = x_24453`.  Both recovered by my own regex over `EQUATIONS.txt`.

**Every one of the 1149 congruence atoms and all 766 off-pin atoms was verified by direct
symbolic expansion through the definition DAG down to `{i1..i6, L, u, P, Q}` and compared to
the formula above.  1149/1149 and 766/766 exact, ZERO mismatches** (`w_verify.py`, `w_final.py`).
This is a recomputation, not a re-reading of P's expansion.

Because `u` is free, the exact integer condition is `c*P | a*L*Z`, i.e. by CRT
**`P | L*Z_k`  AND  `c_k | a_k*L*Z_k`**.  The second is a genuine small-modulus integer
condition invisible mod P; **it only ever restricts, it can never create a family.**

Measured over all 383 blocks: `gcd(a_k,P) = gcd(c_k,P) = 1` for all 1149; **288 congruences
carry `|c|>1` on the handle side and 287 carry `|a|>1` outside** — so P's "the real condition
is `c*P | R`" is confirmed and localised.

## 1.  Nothing escapes through the gate — and booleanity is not even needed

* The liveness cone below the 383 gates is a **pure boolean circuit**: 2,551 nodes —
  256 leaves (**all 256 have an explicit booleanity atom**), 128 constants, 638 aliases,
  765 ANDs, 382 ORs realised as `(a+b) - a*b` on the *same* pair (the 382 SUM nodes occur
  **only** inside those ORs).  Every gate is an AND of two provably-boolean nodes, so
  **L ∈ {0,1} at all 383 blocks** (`w_live3.py`).  My first pass (`w_live2.py`) wrongly
  reported 153 non-boolean gates — that was a **decorrelated abstraction** treating the two
  arms of `(a+b)-a*b` as independent.  Recorded as my own error, caught and fixed.
* **The theorem does not need booleanity.**  For *any* integer `L`:
  * `L ≡ 0 (mod P)`  → congruences vacuous, but the **off-pins force `i5 ≡ i6 ≡ 0 (mod P)`**;
  * `L ≢ 0 (mod P)`  → `P | Z_k` for k=1,2,3, and the 3×2 matrix has rank 2 mod P, so
    `N1 ≡ N2 ≡ 0 (mod P)`; if additionally `L ≢ 1` the off-pins pin the output to 0 as well.
  So **either the law holds mod P or the output is 0 mod P.  There is no third door.**
* Rank 2 mod P is **exhaustive, not sampled**: all 383 matrices are distinct, all six 2×2
  minors per block are nonzero, and `max|minor| = 260,582,651,319,840 < 2^48 << P`, so no minor
  can be a nonzero multiple of P (`w_rank.py`).
* The gate/mux alignment is exact: in **383/383 blocks the output pair `(i5,i6)` is consumed by
  a multiplication gate and that gate is exactly `L`, with no other** (`w_gate.py`).  Every
  short atom touching an output falls into exactly 7 kinds; nothing else reads `i5`,`i6`.

## 2.  THE THEOREM (per gadget, over ℤ with the handles free)

> Write `A = i1-i2`, `B = i4-i3`.  Then modulo P:
> ```
> d(N1,N2)/d(i5,i6) = [[A^2, 0], [B, A]] ,   det = A^3
> ```
> * **`A ≡ 0` forces `B ≡ 0`** (N1 collapses to `-B^2`, and 𝔽_P is a field), and then
>   `N1 ≡ N2 ≡ 0` **identically**: the output `(i5,i6)` is completely free.  — DEGENERACY.
> * **`A ≢ 0`** forces `λ := B/A`, `E ≡ λ^2`, `i3+i6 ≡ λ(i2-i5)`: the output is **uniquely
>   determined** by `(i1..i4)`.  — CHORD.
>
> **There is no third case, and `A ≡ 0, B ≢ 0` is not merely unreachable, it is impossible.**
> "Output has any freedom at all" ⟺ `A ≡ B ≡ 0` ⟺ the gadget sees two equal live inputs.

**Exhaustive machine check of the case analysis** (`w_class.py`): over `p ∈ {5,7,11,13}`,
**every** `Q ∈ 𝔽_p` and **every** `(i1..i6) ∈ 𝔽_p^6` — all `p^7` tuples covered (`N1` does not
involve `i6`, so the loop is `p^6` outer with the `i6` sweep entered exactly when `N1 = 0`) —
the solution set is *exactly* `p^5` degeneracies + `p^4(p-1)` chords, e.g. 371,293 + 342,732
at `p = 13`.
**`THIRD FAMILY` count: 0.  `A=0,B≠0` count: 0.  At every p.**

Over ℤ the classification acquires one extra clause and no extra family: `A ≠ 0` needs
`A^2 | B^2`, which for integers forces `A | B` (verified by enumeration, 0 counterexamples),
so `λ ∈ ℤ`; otherwise the block is **integrally infeasible** — an obstruction, not a solution.
Note the coordinate map `(i1..i6) → (A,B,E,i3+i6,i2-i5)` has all its 5×5 minors in `{0, ±3}`,
so it is only surjective away from characteristic 3; I therefore ran the exhaustive check in
the **original `i`-coordinates**, which sidesteps this entirely.

## 3.  END-TO-END against the deliverable — and a correction to the fleet's account of it

`w_deliv.py`, `w_lie.py`.  The 39,026 witness driven through the full forward map (score
re-derived: 39,026) and every block evaluated **exactly over ℤ**:

| | |
|---|---|
| blocks with all four inputs 0 | **367** |
| gate off, inputs live | **15** |
| **degeneracy `A ≡ B ≡ 0 (mod P)`, gate on** | **1** — block `E = x_33469` |
| chord | 0 |
| **LAW VIOLATED** | **0** |
| nonzero congruence atoms (of 1149) | **0** |
| nonzero leaf pins (of 512) | **0** |

Exactly one degenerate block — **U's §6 finding reproduced from a completely independent
route** (U decoded the curve; I only evaluated `A` and `B`).

> ### REFINEMENT of "it pays 7 equations for a lie on a leaf"
> **It does not lie on a leaf — all 512 leaf pins hold.**  The seven broken atoms are
> ```
> a35759  5113045*(x_7075 * x_9118) - x_29854      OFF-PIN i5 of block E = x_7181
> a35761  x_7075 * x_8731 + x_31864                OFF-PIN i6 of block E = x_7181
> a22229, a22230, a35758, a35760, a35762           the P*u handle/alias atoms of the
>                                                  four corrupted variables
> ```
> All seven sit at the block `E=7181` site: **two off-pins of a DEAD block (gate `L=0`) whose
> output has been driven `≢ 0 mod P`**, plus the five `P*u` handle atoms of the four corrupted
> variables.  The price splits **5 + 2**, not 7 at one place.
>
> **DO NOT read this as "the price of the degeneracy" — I first wrote that and §6a refutes it.**
> Block `E=33469`'s coincidence is **equation-disjoint from all seven** and costs nothing.  The
> 7181 damage and the 33469 degeneracy are two independent things in the same witness.

## 4.  THE LEAD — and it reconciles my own round-1 region

The two off-pins of a block touch a small set of equations.  Over all 383 blocks
(`w_price.py`) the histogram runs 9…20, and the **minimum, 9, is attained by exactly five
blocks**:

| block | off-pin atoms | its 9 equations |
|---|---|---|
| **E=7181** (the deliverable's) | 35759, 35761 | `6816, 8124, 9123, 9421, `**`12231, 12270, 12350, 14584, 18673`** |
| E=3227 | 36124, 36126 | 7081, 11690, 12051, 12233, 17743, 20277, 24141, 32935, 33806 |
| E=4429 | 25138, 25140 | 3068, 7587, 15247, 17392, 24422, 25097, 30542, 31606, 32294 |
| E=30886 | 7516, 7518 | 658, 3005, 4489, 13891, 15141, 15635, 17675, 30993, 33618 |
| E=31606 | 31199, 31201 | 4655, 7223, 12086, 16607, 17668, 18924, 23660, 29322, 35517 |

**All ten pairwise overlaps are 0.**  The four non-deliverable ones share nothing with the
failing set.

**This derives my round-1 frame-B region from scratch.**  Block 7181's nine equations are
five of the seven failures **plus `6816, 8124, 9123, 9421` — four of the six essential rows
I found last round** (`{2554, 6816, 8124, 9123, 9421, S}`), which are also four of the six
"prices" in the 32-way trade.  **K (34 free inputs) IS the off-pin neighbourhood of block
7181**, arrived at here by pure structure with no linear algebra.

## 5.  SCOPE — stated as the rules require

| claim | status |
|---|---|
| the 5-atom integer form of every block | **exhaustive**, 1149+766 symbolic identities, 0 mismatches |
| rank 2 mod P of the 3×2 matrix | **exhaustive**, all 383, via minor magnitudes |
| gate ∈ {0,1}; gate/mux alignment; off-pins present | **exhaustive**, 383/383 and 766/766 |
| the two-family classification | **proved** over any integral domain; the *case analysis* additionally **machine-checked exhaustively** over 𝔽_p, p ∈ {5,7,11,13}, all Q, all (i1..i6) |
| "no third family" | **holds at the ATOM level.**  See the boundary below. |
| the off-pin incidence table (9…20) | **exhaustive** over 383 blocks — but it is an *incidence*, **a screen, not a price**, and §6b shows it does **not** predict price |
| injection cost at the 5 minimum-incidence blocks | **budget, not exhaustion** — one magnitude, 1–2 slots per block (§6b) |
| `minbreak(P) = |P|`, gain 0, on the corrected knob set `K+ = 40` | **exhaustive over the essential-row family**; general breaks are **not** re-proved on `K+` (§6c) |

> ### THE ONE BOUNDARY, stated plainly
> # **THE CLASSIFICATION IS CLOSED AT ATOM LEVEL AND OPEN AT EQUATION LEVEL.**
> Everything in §2 classifies the solutions of **`atoms = 0`**.  It does **not** classify the
> solutions of `equations = 0`, which is a strictly larger set.  Do not cite §2 as "the gadget
> has only two solution families" without this sentence attached.
>
> The checker requires each **equation** to vanish, and an equation is a coefficiented sum of
> ~12 atoms (congruence atoms sit in **9–16 equations each, mean 12.28; none is ever alone in
> an equation**).  Everything above classifies the solutions of **atoms = 0**.  Equation-level
> cancellation between atoms is a strictly larger solution set and this theorem does not cover
> it — that is exactly the trade machinery I measured in round 1 (32-way, 1-for-1, gain 0
> inside K).  **The deliverable itself does not use it at gadget level: all 1149 congruence
> atoms and 764 of 766 off-pins evaluate to exactly 0 in the witness.**

## 6.  THE OUT-OF-K MEASUREMENT — done, and it is a NEGATIVE with two corrections

### 6a.  CORRECTION to my own §3: the degeneracy is FREE; the 7 are a SEPARATE obstruction

`w_chain.py`, `w_chain2.py`.  Block `E=33469`'s four live inputs are themselves **free inputs**
in this orientation, set directly to coinciding values.  The check atoms they move touch **46
equations, and the overlap with the seven failures is EMPTY** — as is the overlap with the two
broken off-pins.

> **The degeneracy at 33469 costs nothing.**  The seven failures live entirely at the block
> `E=7181` site.  My §3 accounting (2 off-pins + 5 handle atoms) is right; the *causal* reading
> "it pays 7 equations for the degeneracy" — mine as much as the fleet's — is **wrong**.  They
> are two independent things that happen to sit in the same witness.

Scope: this is a statement about the witness and its 40-knob neighbourhood, **not** a claim
that the 7181 site is unsatisfiable in principle.

### 6b.  The four equation-disjoint blocks are all MORE expensive, not less (`w_sites.py`, `w_outk.py`)

All five minimum-incidence blocks have a **dead gate** (`L = 0`); at the four non-deliverable
ones every slot is `≡ 0 mod P`, so their off-pins currently hold.  Injecting (driving a free
output slot to a value `≢ 0 mod P`), priced **exactly through `frameB.State`**:

| block | free output slots | NEW equations broken | score |
|---|---|---|---|
| **E=7181** (the deliverable's) | 9118, 8731 | **5** (of its 9) | 39,026 |
| E=30886 | 18957, 6120 | **8** | 39,018 |
| E=3227 | 36247, 26738 | **9** | 39,017 |
| E=4429 | 11131, 35676 | **9** | 39,017 |
| E=31606 | 15317, 9121 | **10** | 39,016 |

> **The 9-equation incidence floor is NOT the price.**  Across the five blocks that share it
> the actual cost runs 5…10.  The ranking that produced them is a **screen, and a weak one**;
> the deliverable's site is the cheapest of the five by three equations.  **No score movement.**
>
> **Budget, not exhaustion:** one injection magnitude and 1- or 2-slot combinations per block.
> A cleverer value or a wider slot set could be cheaper; I did not enumerate.

Control, same file: zeroing 7181's own outputs gives **39,023 / 39,022**, i.e. strictly worse —
the injected values are load-bearing for `{6816, 8124, 9123, 9421}`.  That is the k-for-k trade
of round 1 seen directly, without the linear model.

### 6c.  A GAP IN MY OWN ROUND-1 SCOPE, found and closed (`w_setup2.py`, `w_exhaust2.py`)

Round-1's `K` was *"free inputs reaching a **nonzero** atom, ∪ carriers of S"*.  But an equation
is a sum of atoms: to repair a broken one you may move a **different, currently-zero** atom in
it to compensate — and those movers were baked in as **constants**.  **So K was not the
complete mover set, and I did not notice.**

Measured: the 7 failing equations contain 24 atoms; their movers number **40**, so round-1's 34
missed **6** — `[4068, 11368, 16763, 21574, 24559, 32230]`.  Re-running the setup on `K+ = 40`:
205 rows, 198 satisfied, `rank(A_SAT) = 32` (was 26), and the essential rows are **the same six**
`{2554, 6816, 8124, 9123, 9421, S}`.

> Re-ran round-1's exhaustive test on the extended system — all 2⁶ essential break-sets × all
> 127 bought-sets, exact integer oracle, 105 s:
> **`minbreak(P) = |P|` exactly for every `|P| ≤ 6`, GAIN = 0 everywhere, all seven unbuyable.**
> **The six missing movers add nothing.**  The omission was real and it was inert.

**One thing does NOT carry over:** on `K+` the packing lemma yields only `t = 1` disjoint
full-rank subsets of the 192 redundant rows (round 1 had enough for `b ≤ 2`).  So the `K+`
result is exhaustive **over the essential-row family only** — the general-break exhaustiveness
of round 1 (cocircuits, `b ≤ 2`) is a statement about `K = 34` and has **not** been re-proved on
`K+ = 40`.

### 6d.  Where I would go next

1. **Redo the cocircuit / general-break closure on `K+ = 40`**, so the `b ≤ 2` exhaustiveness
   is restored on the corrected knob set.  This is now the top item: round 1's strongest
   exhaustive row is currently scoped to a knob set I have shown to be incomplete.
2. The `s = 3..6` cocircuit gap (round-1 item #1) — unchanged, still open.
3. **Not** more injection sites.  The five-block screen came back negative and the incidence
   ranking that generated it does not predict price; a broader sweep would need a real price
   oracle per block, which is items 1–2's machinery anyway.

## 7.  ITEM 1 — the cocircuit closure redone on `K+ = 40`.  **IT DOES NOT CLOSE AT ROUND-1
## STRENGTH: round 1's `s = 3..6` search was materially incomplete, and I can now say by how much.**

### 7a.  Preconditions re-verified on `K+` first (`w_hom2.py`)
All 198 SAT rows are **homogeneous** (rhs = 0 for every one), `rank_Q(A_SAT) = 32`, and the full
205-row system is still **ℚ-consistent** (`rank = rank[A|b] = 34`).  So round 1's structural
result — *the whole obstruction is integrality, only rank-dropping deletions matter, and any
break-set is a union of minimal cocircuits* — **reproduces verbatim on the corrected knob set.**

### 7b.  A new enumerator with NO window and NO skipped cases (`w_cocirc3.py`)
Round 1 (and my first port, `w_cocirc2.py`) chose `s-1` outside columns from an adaptive window
and *skipped* the subsets where they were dependent — 3.07M skips in round 1.  The replacement
is an exact case split:

> a codeword with `|supp| ≤ 6` is nonzero on at most `6-s` outside columns, so **among ANY
> `7-s` outside columns at least one vanishes at it.**

Branch on that: pick `7-s` columns whose functional is not identically zero on the current
subspace `V`, recurse into `V ∩ ker φ_k` for each (every cut drops the dimension by exactly 1).
Depth `s-1`, branching `7-s`, so the tree is finite and **complete**.  When fewer than `7-s`
columns are live on `V`, *every* point of `V` already has weight `≤` that count — the whole
subspace is **recorded** instead of skipped, which is exactly the case round 1 threw away.

### 7c.  How much round 1's window was missing — measured head to head

| | window (`w_cocirc2.py`, same K+) | **exact (`w_cocirc3.py`)** |
|---|---|---|
| supports at `s ≤ 3` | 1,485 | **3,184** |
| supports at `s ≤ 4` | 5,965 | **25,473** at `s=4` alone |
| positive-dimensional nodes | **skipped** | **recorded**: 2,376 at `s≤3`, 65,619 at `s=4` |

> **The window method found well under a quarter of the `s ≤ 4` supports.**  Round 1's
> `j = 4..7` row was labelled *budget, not exhaustion* and that label was right — but the gap
> was **much larger than the 3.07M skip count suggested**, and it is not closed by re-running
> the same method on more information sets.

### 7d.  What the exact enumeration gives, so far (`w_close3_s13.out`)
From the `s ≤ 3` exact data alone: **60 minimal cocircuits**, sizes `{1:6, 2:2, 3:2, 4:3, 5:7,
6:40}`, **every one re-verified genuinely rank-dropping over ℚ (0 mod-p artefacts)**; union
closure to size ≤ 6 gives **510 break-sets**.

* The six size-1 cocircuits are the **same six essential rows** `{2554, 6816, 8124, 9123, 9421, S}`.
* There are **two** minimal size-2 cocircuits — `{22563, 8687}` (round 1's) **and `{36489, 8985}`,
  which round 1 on `K = 34` never saw.**

### 7e.  Status and what is still running (PIDs recorded)
| job | file | PID | state |
|---|---|---|---|
| exact `s = 1..3` | `w_cocirc3_raw_s1_3.json` | — | **done** |
| exact `s = 4..6` | `w_cocirc3_raw_s4_6.json` | 28998 | running (`s=4` printed at 277 s) |
| exact `s = 4` only | `w_cocirc3_raw_s4_4.json` | 29354 | running (persists `s=4` without waiting for `s=6`) |
| window, 3 information sets | `w_cocirc2_raw.json` | 9543 | running (set 0 through `s=5`) |
| closing test on exact `s≤3` | `w_close3_s13.json` | 31769 | running |

**The box is at load ≈ 21 on 4 CPUs** (other agents), so these are far slower than their solo
timings (`s=4` exact: 277 s of its own time).

> **Honest status of the frame-B optimality row:** `minbreak(P) = |P|`, gain 0 is confirmed on
> `K+` **over the essential-row family** (§6c, exhaustive) and over the **510** break-sets from
> the exact `s ≤ 3` enumeration.  **General breaks at `j = 4..7` remain BUDGET, not exhaustion**,
> and after §7c that budget should be described as *substantially* short, not marginally.
> Finishing `s = 4,5,6` exactly is what converts it — the enumerator now exists and is correct;
> it only needs CPU.

## Re-entry
```
cd solve_lab/agentW_work
python3 model.py ; python3 fwd2.py     # rebuild the pkls (~30 s)
python3 w_blocks.py ; python3 w_cong.py ; python3 w_ring.py   # the 383 blocks + the ring
python3 w_rank.py      # rank 2 mod P, exhaustive                      (~1 s)
python3 w_live3.py     # the gate is boolean                           (~30 s)
python3 w_gate.py ; python3 w_offpin.py                                (~1 min)
python3 w_verify.py    # 1149 congruence identities, symbolic          (~52 s)
python3 w_final.py     # 766 off-pin identities + handle privacy       (~15 s)
python3 w_class.py     # THE CLASSIFICATION + exhaustive small-field   (~2 min)
python3 w_deliv.py ; python3 w_lie.py  # end-to-end on the 39,026 witness
python3 w_price.py     # the five minimum-incidence blocks
python3 w_chain.py ; python3 w_chain2.py   # the degeneracy is equation-disjoint: costs 0
python3 w_sites.py ; python3 w_outk.py     # the four out-of-K blocks; K vs K+
python3 w_exhaust2.py  # round-1's exhaustive test re-run on K+ = 40   (~105 s)
python3 w_hom2.py      # K+ preconditions: homogeneous, rank 32, Q-consistent
WSMIN=1 WSMAX=3 python3 w_cocirc3.py   # EXACT cocircuits, no window, no skips
WFILES=w_cocirc3_raw_s1_3.json WOUT=w_close3_s13.json python3 w_close3.py
```
Artifacts: `w_blocks*.json`, `w_verify.json`, `w_class.json`, `w_deliv.json`, `w_price.json`.

---
---

# ROUND 1 (check-in 92) — preserved verbatim below

**Best verified score: 39,026 / 39,033** — same file, same 7 failures.  **I did not beat it.**

## TASK 2 — SETTLED.  O's `|K| = 34` is CORRECT in frame B.
Frame B = `frameB.Frame([642, 28730, 29854, 31864])` reproduces the witness bit-for-bit:
score **39,026**, same 7 failures, **0 of 38,748 variables differing**, 7 nonzero check atoms
`[22229, 22230, 35758, 35759, 35760, 35761, 35762]`.

| orientation | free inputs | witness score | nonzero atoms | \|U\| | \|C\| | \|U∩C\| | **\|K\|** |
|---|---|---|---|---|---|---|---|
| **frame B** `[642,28730,29854,31864]` | 8,751 | 39,026 | 7 | **15** | **26** | **7** | **34** |
| default (no detach) | 8,747 | **39,020** | 5 (incl. a37887) | 30 | 26 | 26 | 30 |

T's rebuild (12 / 11 / 23, overlap 0) is a fact about **F's parse in the default orientation**.
Ledger row → CONDITIONAL, scope: agent H's model, frame B's orientation.  **Not a defect.**

## TASK 1 — RESULTS (round 1)
* **(a)** the 21 triples O never reached: all `b<=2 exhausted, none`; 298,158 exact integer
  solves, 2,065 s (`w_j3.py`, `w_j3.log`).  O's brief said 14 = "every triple containing
  eq12231"; there are 15, so the gap was **21**, not 20.
* **(b)** the linear model's pricing is EXACT outside the model: `w_trade_12231_break2554.json`
  → checker 39,026/39,033, 27 variables from the witness.
* **(c)** AUDIT: **the trade is 32-way, not 7-way.**  eq8680 is one of six prices and the
  unique price for **eq29125 alone**.  O's collateral accounting is sound; O's Lemma untouched.
* **(d)** `rank([A|b]) = 28` over all 175 rows with the rhs **not** a pivot → **the full system
  including all seven failing rows is CONSISTENT over ℚ.  The entire frame-B obstruction is
  integrality; no ℚ or LP relaxation can ever prune here.**  Exactly **6 essential rows**:
  `{2554, 6816, 8124, 9123, 9421, S}`.
* **(e)** exhaustive over the essential-break family, all j=1..7: **`minbreak(P) = |P|`
  exactly, gain 0 everywhere**; all seven unbuyable at any `b <= 6`.
* **(f)** honest correction to my own (e): redundant-row breaks are **not** worthless —
  `{22563, 8687}` is a genuine minimal cocircuit containing no essential row.
* **(g)** cocircuits: 70 minimal, sizes `{1:6, 2:1, 3:2, 4:5, 5:14, 6:42}`.  `s<=2` **exact**;
  `s=3..6` is a **bounded search** (3,070,206 degenerate subsets skipped).
* **(h)** 520 union-of-minimal-cocircuit break-sets × 127 bought-sets, 6,806 exact solves,
  49 s: **BEST GAIN = 0.**

### Round-1 scope
| claim | status |
|---|---|
| j=1 b=0 ; j=2 b<=1 ; **j=3 b<=2 all 35 triples** | **exhaustive** (the last twice: brute force and structurally) |
| j=1..7 restricted to the 6 essential rows | **exhaustive** |
| j=4..7 general breaks | **budget, not exhaustion** |
| all of it | scope: **34 of 8,751 free inputs**, frame B, agent H's model |

### Do not redo
* The ℚ relaxation of the frame-B system — vacuously YES, nothing can prune.
* Brute-forcing break-sets that do not drop `rank_Q` (162 of 168 rows are redundant).
* Reading a greedy net-zero as a negative.
