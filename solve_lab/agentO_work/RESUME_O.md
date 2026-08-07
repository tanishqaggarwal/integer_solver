# RESUME_O — agent O.  Consolidated and closed.

**Best score: 39,026 / 39,033.  I did not beat it.**  Three artifacts of mine verify at that
score with `solve_lab/checker.py` (all my values are ≤ ~3,050 bits, so plain `checker.py`
parses them; `verifyE.py` was never needed):
`grow23618_39026.json`, `region_opt_39026.json` (a *distinct* point — different values for all
seven region variables), and the witness itself, re-verified:
failing `[12231, 12270, 12350, 14584, 18673, 22044, 29125]`.

Read §1 first; it is the one result other threads depend on.

---

## §1. THE LEMMA — `S = 0` is forced.  Unconditional, audited, citable.

> **`eq8680`'s left-hand side equals `S⁴`, where `S` is an AFFINE form in atoms.**
> `checker.py` requires the LHS to be exactly `0` over **ℤ**, an integral domain, so
> `S⁴ = 0 ⟺ S = 0`.
> **Therefore `S = 0` in every satisfying assignment.**
> No knob set, no frame, no configuration, no divisibility, no modulus.

**Cite this section; do not re-derive it.**  It is load-bearing in three other threads: it is
why δ₀ died (§4), it is the mechanism behind N's finding that eq8680 is exactly what detaching
`x_28730` buys, and it is why M's enumeration space shrank from 2¹⁸ to 2¹⁶.  Per N it is also
precisely the 39,025 → 39,026 step (N's result, not my measurement).

### How it was established, and how it was audited
First stated from H's parse, where eq8680 has *one* term (`a37887`, source literally `(X)·(X)`).
Agent T audited it and found three errors.  I re-verified all three **against the raw text of
`EQUATIONS.txt`, using no parser at all** (`verify_lemma.py`, `runs/verify_lemma.log`) — the only
way to settle a disagreement between two parses.  Perturb one variable, read the raw LHS:

| perturbation | `S` | raw LHS | `S²`? | `S³`? | `S⁴`? |
|---|---|---|---|---|---|
| `x_4432 += 2` | 2 | 16 | ✗ | ✗ | ✓ |
| `x_4432 += 3` | 3 | 81 | ✗ | ✗ | ✓ |
| `x_4432 += 5` | 5 | 625 | ✗ | ✗ | ✓ |
| `x_19964 += 2` | −2 | 16 | ✗ | ✗ | ✓ |
| `x_28730 += 3` | −3 | 81 | ✗ | ✗ | ✓ |
| `x_23754 += 2` | −18 | 104976 | ✗ | ✗ | ✓ |

`LHS == S^k` for **k = 4 only**.  Corrections, all confirmed:
1. **`S⁴`, not `S²`** — the nesting is two levels, `LHS = T·T` with `T = S·S`.
2. **The object with slope +1 is `S`**, the affine form: `dS/dx_4432 = +1`,
   `dS/dx_19964 = −1`, `dS/dx_28730 = −1`.  I had written "`eq8680 = T²`, `T` linear" *and*
   quoted `dT/dx_4432 = +1`; those cannot be the same object (`T = S²` ⇒ `dT/dx = 2S+1`).
3. **18 vs 20 is granularity, not contradiction** — the raw text has **18** bracketed groups;
   E emits **20** `(coef, atom)` entries because it splits exactly two of them,
   `−13·(x_21279·x_31731 + x_35619)` → a23622,a23623 and
   `−5·(x_34600 − x_30108 + x_23642)` → a11876,a11877.  **18 + 2 = 20**, both correct.
   ⚠ **`S`'s 18 is NOT M's enumeration exponent 18.  Different 18s.**

**The error was in the prose alone.**  My frame-B "S row" was built from H's *inner* factor,
which **is** this affine form — I measured its slope as +1 before using it — so every search
below constrained the right object.  *This number is wrong* ≠ *this result is wrong*.

The conclusion is robust to the exponent entirely: `S^k = 0 ⟺ S = 0` for any k ≥ 1.  In **E's**
decomposition `eqfails` tests the affine sum itself, so E's model states `S = 0` directly while
H's states `S⁴ = 0` — same zero locus.

**Cross-link:** all three p-handles T found missing from L's census are terms of `S`, confirmed
by source match — `25·(x_18253 − x_4339·x_15120)` = a20450,
`1·(x_37720 − x_14466·x_35531)` = a20452, `23·(x_23642 − x_8173·x_10422)` = a11875.

---

## §2. The seven-way 1-for-1 trade

**Knob set (inline, so the claim travels with its scope):** frame B,
`frameB.Frame([642, 28730, 29854, 31864])`, which reproduces the witness bit-for-bit
(39,026, same 7 failures, **0 variables differing**).  Twelve knobs —
`{642, 1329, 9413, 10903, 17325, 28730, 29854, 31864}` (region-private)
`+ {7068, 4432, 8731, 9118}` (carriers) — all free inputs there.  They reach exactly
**12 checks / 29 equations**, and **all 7 witness failures are inside**; nothing unreachable.

> **Every one of the 7 failing equations is individually buyable, and every purchase costs
> exactly `eq8680`.**  Score pinned at 39,026 seven ways; the failing set merely rotates.
> No pair is buyable.

**Why:** `a23618 = x_4432 − x_19964 − x_28730` is `S`'s first term at coefficient +1, and
`dS/dx_4432 = +1`, `dS/dx_28730 = −1`, zero for every other region knob.  So `S = 0` is exactly
`δx_4432 = δx_28730`.  With `S = 0` imposed as an explicit row, **nothing is buyable at all**.

This confirms and extends agent H's "a22231 buys 1 row and costs eq8680, exactly" to **all
seven**, and supplies the derivative behind agent G's characterisation of eq8680 as a binary
quadratic form of discriminant exactly 0.

**Scoped theorem.**  Let `U` = the 15 frame-B free inputs reaching any of the witness's nonzero
check atoms, `C` = the 26 carriers of `S`; `K = U ∪ C`, **|K| = 34**.
**Every assignment agreeing with the witness outside `K` satisfies at most 39,026 equations.**
Over `K`: 64 reachable checks, 190 reachable equations, all 7 failures reachable, 175 rows.

**Model exactness was verified before trusting a negative from a linear model:** a 5-point probe
(t = 1,2,3,5,7) finds precisely the same 7 non-affine checks the 2-point probe found, **none
missed**.  The 16 dropped rows all contain one of those 7 and **none currently fails**, so
dropping them is *permissive* — the solver was free to break them and still found nothing,
making the negative strictly stronger.

---

## §3. The compensation channel, and its budget stated as a budget

`S = 0` can be evaded only if some of `S`'s **other terms** move to keep it zero.  That is the
only route by which the 1-for-1 trade could be leveraged.

**There is no free compensator.**  All 20 atoms of `S` live in **10–18 equations** — none is
confined to eq8680 the way `a37887` is in H's bundled parse.  Nine are checks in E's frame.  The
equations they disturb are **the region's own** (2554, 6816, 8124, 9421, 12231, 12270, 12350,
14584, 22044, 29125), so every carrier of an `S` component is already a carrier of `a37887`, and
**all 26 were in `K` from the start.**  *The channel was never a missing knob; it is purely a
budget.*  Full atom table in `T_COMPENSATION.md`.

To beat 39,026 with `j` bought and `b` broken we need `b < j`:

| budget | scope actually tested | solves | result |
|---|---|---|---|
| `j=1, b=0` | **complete** | 7 | none |
| `j=2, b=0` | **complete** (all 21 pairs) | 21 | none |
| `j=2, b≤1` | **complete** — 21 pairs × each of 168 satisfied rows *and* the `S=0` row | 3,570 (21 s) | **none** |
| `j=3, b≤2` | **14 of 35 triples enumerated completely** (per triple: b=0, all 168 b=1, all C(168,2)=14,028 b=2) | 198,772 (33 min) | **none** |
| `j≥4` | greedy upper bound only | — | drops 25–26 vs needing <4 |

All 21 pairs were individually feasible, so nothing was vacuously pruned at `j=2`.

**The trap I nearly fell into.**  The greedy pass flagged `[12231, 12270, 12350]` as dropping
*exactly 3* — net zero.  Greedy only **upper-bounds** the drops, so the true minimum could have
been 2, **i.e. 39,027**.  Enumerated properly at `b ≤ 2`: none.  Reading the greedy number as a
negative would have made this a false negative.

**Stated as budget, not exhaustion:** the 14 completed triples are exactly those containing
eq12231; the other 21 were **not reached** within the wall-clock cap; `j≥4` is greedy only.
So: the 1-for-1 trade is proven unleverageable **at budget 2 over `K`**, and **at budget 3 for
every triple containing eq12231**.  Scope throughout: **34 of 8,751 free inputs, frame B's
orientation**.

---

## §4. The δ₀ line, and why it died

The witness's residual is 8 atoms touched by exactly 12 equations (13 with `a23618`, the extra
one being **eq8680**), with 7 variables private to the region — moving them cannot break
anything else, so **failing equations = |E(R)| − maxsat(R)** exactly.  Over ℚ the system has a
unique solution satisfying all 13; over ℤ **five** coordinates are blocked, by moduli
`p, p, p, 2458959, 2458959·p`.

Rather than sample, I put the boundary shift into the unknown vector and solved `A z + B δ = b₀`
over ℤ.  **0 of 9** single supports work, **0 of 36** pairs, **0 of 84** triples; **12 of 126**
quadruples do — including exactly `{a23616, a23618, a36660, a36662}`, the four constants that
are *not* p-multiples.  Applying δ₀ makes **all 13 region equations hold**, verified end-to-end.
Two of the four carriers are free (`x_8731`, `x_9118` — H's zero-collateral knobs).

**It is not realisable.**  `a23618` is the sole carrier of the `L` shift, and `S = 0` forces
`δx_4432 = δx_28730` — **collapsing the shift direction onto the private handle direction and
annihilating precisely the degree of freedom δ₀ needs.**  My atom-level model could not see this:
it drops `a37887` as nonlinear, and `a37887` is the one atom *outside* the region tying the two
carriers together.  I emitted the handoff before pricing it and then issued
**`DELTA0_STATUS.md`** correcting it, so no one spent cycles on a dead target.
δ₀ remains a correct lattice target *for the region in isolation*; the coupling to `a37887`
kills it.  Exact values in `DELTA0_FOR_M.json`; read `DELTA0_STATUS.md` first.

---

## §5. The rate — computed before spending the cores, and the scan not run

I proposed a 2,800-configuration scan, then computed its cost first.

- **Rate ≈ 2⁻⁷⁶⁷.**  Admissible boundary changes form a coset `δ₀ + Λ₀`; measured period in each
  direction with an exact solvability oracle: **p** for `a23618`, `a36660`, `a36662`; for
  `a23616` the period exceeds every modulus tested (up to 2458959·p).  So `[Z⁴ : Λ₀] ≥ 2^768`.
  Expected hits in 2,800 configurations ≈ **10⁻²²⁷**.
- **Worse than a bad rate — zero variance.**  Across **35 configurations** the four boundary
  quantities `K1, L, K2, J` are **identically 0** — 1 distinct value out of 35, each.  The scan
  would have measured a single point 2,800 times.  They are **assignment knobs, not
  configuration knobs**: the witness has `J ≠ 0` because it *assigns* free variables.

**Do not run that scan.**

---

## §6. Earlier rounds, condensed
- **Complete singles sweep**, all 106 proven pin-solvable bits at E's cfg0 (`runs/singles.jsonl`):
  none beats the empty set.  "Only 2 representatives per channel" was not the explanation.
- **The monotonicity tension, resolved by mechanism.**  cfg0 is the **(0,1) branch**
  (`x_7715 = 0`, `x_34554 = 1`, so `x_15298 = 0`), which makes atoms 20649/20652/32148 *vacuous*.
  Turning on **any** a-tree bit sets `x_15298 = 1` and fires all three at once.  One gate is the
  whole monotone cost — not per-bit accumulation.
- **E's frame provably cannot represent the witness**: feeding the witness's own free values
  through `engine.forward` gives 25 failing equations, not 7, with four divergence roots
  `x_31864, x_29854, x_642, x_28730` — the same four M and P identified.
- The witness's repair is **outside the dependency cone of the residual**: cone(20649,20652,
  32148) has 277 free variables and **none** of the six carriers; each carrier's +1 probe has
  zero delta on all three rows.  No cone-generated closure at any `maxr`/`maxv` can reach it.

---

## §7. Re-entry and files
```
cd solve_lab/agentO_work
python3 verify_lemma.py     # re-verify §1 against raw EQUATIONS.txt, no parser  (~1 min)
python3 fb_max.py           # the seven-way 1-for-1 trade                         (~2 min)
python3 fb_sq.py            # with S = 0 imposed: nothing buyable                 (~2 min)
python3 fb_U.py             # scoped theorem over K, |K| = 34                     (~3 min)
```
Key documents: **`EQ8680_LEMMA.md`** (§1, corrected + audited), **`T_COMPENSATION.md`** (§3),
**`DELTA0_STATUS.md`** (§4 — read before `DELTA0_FOR_M.json`).  All logs in `runs/`.
`agentH_work` was imported **read-only** throughout; 0 files modified there.  No git commands.

## §8. Do not redo
- The 2,800-configuration scan (§5) — dead twice over.
- Singles over all 106 representatives at cfg0 — done, negative.
- Deeper closure on E's residual (`maxr` 3→6, `maxv` 2000→8000, unfiltered knobs, lazy rows):
  never below 28 failing.  The knobs were never the problem — the **cone** is.
- Single- and double-atom growth of the witness region: all fail on eq29125 mod p.
- `(a,b)` pairs through `full11.solve_pair`: the first 68 all give **39,013** with residual
  exactly `{20649, 20652, 32148}`, independent of `b`.

## §9. What is left, honestly
Within `K` the region is closed to the budgets in §3.  The untested remainder: the **21 triples
not containing eq12231** at `j=3, b≤2`, and `j≥4` beyond greedy.  Both are mechanical extensions
of `fb_j3.py` (raise `OCAP`).  Beyond that, any improvement must come from **outside `K`** —
i.e. from an assignment that differs from the witness in one of the other 8,717 free inputs,
which §1 does not constrain.
