# eq8680 — the constraint is forced, and it is LINEAR.  Answer: option (2).

Agent O.  Coordinator asked which of two it is: can `S ≠ 0` be reached at acceptable cost, or is
`S = 0` genuinely forced?  **It is forced**, and the argument is parser-independent.

## The Lemma, stated correctly

My first version of this was wrong in a way worth recording.  In agent H's parse, eq8680 has
**exactly one term** — coefficient 1 on atom `a37887`, whose source is literally `(S)·(S)` — and
I was about to conclude "there is no other atom to compensate with".  That is an artifact of how
H's parser *bundles* the expression.  **E's independent parser sees the same equation as 20
terms with `issq = True`.**  Had I claimed the theorem from H's parse alone it would have rested
on a bundling choice.

Both parses agree on the mathematics, and the correct statement is:

> **Lemma.** `eq8680 = T²`, where `T` is a **linear form in atoms**:
>
> `T = a23618 + 6·a23619 + 15·a23620 − 21·a23621 − 13·a23622 − 13·a23623 + 25·a20448 + a20449
>    + 25·a20450 + 28·a20451 + a20452 − 4·a20453 + 23·a11875 − 5·a11876 − 5·a11877 + 20·a11878
>    − 27·a11879 + 35·a11880 + 17·a11881 − 14·a11882`   (E's atom numbering)
>
> Verified in two independent parses: E's (20 atoms, `issq = True`, outer multiplier 1) and H's
> (one atom `a37887 = T²`, `sq = True`, multiplier 1, and the two factors of the source are
> byte-identical).  In H's frame `T = S`.
>
> **Therefore every assignment satisfying eq8680 has `T = 0`.**

The reason nothing can compensate is not that the equation has one atom — it is that **a square
has a single zero locus**.  `eq8680 = T²` vanishes only on the hyperplane `T = 0`.  This is
unconditional: no knob set, no frame, no divisibility.

## Why this is exactly the obstruction to δ₀

`a23618 = x_4432 − x_19964 − x_28730` enters `T` with coefficient **exactly +1**, and it is the
sole carrier of the `L` boundary shift δ₀ requires.  Measured (exactly, at t = 1,2,3,5,7):
`dT/dx_4432 = +1`, `dT/dx_28730 = −1`, zero for every other region knob.  So `T = 0` forces
`δx_4432 = δx_28730` **unless one of the other 19 atoms in `T` moves to compensate** — and those
19 atoms are the precise, explicitly-named compensation channel.

## The scoped optimality theorem

> Let `U` be the frame-B free inputs reaching any of the witness's nonzero check atoms
> (**|U| = 15**), and `C` the carriers of `T` (**|C| = 26**); `K = U ∪ C`, **|K| = 34**.
> **Every assignment agreeing with the 39,026 witness outside `K` satisfies at most 39,026
> equations.**

Evidence, all exact:
- The 7 failing equations depend only on the witness's nonzero atoms — verified term by term.
- Over `K`: **64** reachable checks, **190** reachable equations, **all 7 failures reachable**,
  none permanently failing.
- **175 rows, every one exactly affine.**  Validated by a 5-point probe (t = 1,2,3,5,7): it
  finds precisely the same 7 non-affine checks the 2-point probe found, **none missed**, so the
  linear model is exact and the negative result is sound.  The 16 dropped rows all contain one
  of those 7 checks and **none of them currently fails**, so dropping them is *permissive* —
  the solver was free to break them and still found nothing, making the negative strictly
  stronger.
- **Zero-collateral:** no subset of the 7, of any size, buyable while keeping all 168 satisfied.
- **Net gain, pay 1** (drop the `T = 0` row, i.e. let eq8680 break): no subset of size ≥ 2.
- **Net gain, pay 2** (eq8680 + any one other satisfied row): no subset of size ≥ 3 or 4.

So the trade is exactly 1-for-1 and cannot be leveraged: **pay one, buy at most one.**

## What this is NOT
- It is **not** a global optimality proof.  The scope is "agrees with the witness outside 34 of
  the 8,751 frame-B free inputs".  An assignment differing elsewhere is untouched by it.
- Test B was budget-capped (1500 s); it covered k = 3 and 4 against every single extra payment,
  not every pay-2/pay-3 combination.
- It assumes frame B's orientation, which reproduces the witness bit-for-bit (0 variables
  differing) but is one orientation among several.

## Where the remaining freedom is, precisely
The 19 atoms of `T` other than `a23618`, and their carriers.  A move that changes `T`'s other
19 atoms so as to hold `T = 0` while `δx_4432 ≠ δx_28730` restores the `L` direction δ₀ needs.
All 26 carriers of `T` were in the knob set above and none sufficed **singly or in the
combinations the solver explored** — but the solver was constrained to keep the other 168 rows,
so the open question is whether a *deliberately budgeted* multi-atom compensation inside `T`
exists.  That is the one door I did not close.
