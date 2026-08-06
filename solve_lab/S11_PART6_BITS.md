# Session 11, Part XVII — the obstruction as conserved quantities of the message

Part XIV reduced the instance to a GF(p) circuit and showed no continuous move clears it.  This
part turns that negative into a usable positive: the obstruction certificates are **invariants**,
computable from any message in 0.08 s, and they factor over the OR-trees.

## 1. Certificates are conserved quantities

A certificate y satisfies `y . J = 0` over every continuous knob, so

        INV_y  =  sum_a  y_a * r_a          (r = check residues)

is **constant under every continuous move**.  Verified directly: 39 of 40 random single-knob
perturbations leave all six invariants bit-for-bit unchanged (`s11/bits10.py`).

A full solve needs every r_a = 0, hence needs every invariant to vanish.  So the search objective
stops being "how many checks fail" — a signal that plateaued at 4 and gave no gradient — and
becomes **"drive six explicit GF(p) numbers to zero"**, each evaluable in one forward pass.

## 2. They factor over the OR-trees

Measuring which bits move which invariant, from three different anchors (`s11/bits15.py`):

    anchor {x24601}  (U=1,V=0)   inv0,1,4: 64 bits (D 41, C 23)   inv5: 18 bits (all C)
                                 inv3: 166 (A,C,D)   inv2: 241
    anchor {x2081}   (U=0,V=1)   inv0..4: 255 bits                inv5: 36 bits (all C)
    anchor {x2081,x24601}        inv0..4: 256 bits                inv5: 37 bits (all C)

**Invariant 5 is moved only by C-tree bits, at every anchor** — it is a function of the C-subset
alone, independent of the other 219 bits.  And in the channel-C regime (V=0), invariants 0, 1 and
4 depend only on the V-side (C u D, 78 bits).  That splits a 2^256 problem into staged pieces.

## 3. Invariant 5 is highly degenerate

Enumerating C-subsets of weight <= 3 (`s11/bits16.py`), inv5 takes only ~300 distinct values over
3,000 subsets — nothing like a random GF(p) function.  Whatever it measures, it collapses.

## 4. The bit landscape, mapped exactly

`s11/bits5.py`, `bits7.py`, `bits8.py` score messages exactly rather than perturbing one pattern:

* **All 256 weight-1 messages**: 15-30 failing checks.  Best is `{x24601}` at 15 (channel C),
  then `{x2081}`, `{x4287}`, `{x13195}` at 17 (channel B).
* **Exhaustive weight-2** (all 32,640): the minimum is **4**, reached by exactly three messages —
  `{x24601, x2081}` (the checkpoint), `{x24601, x4287}`, `{x24601, x13195}`.  All three fail the
  *same four atoms* a7930, a29539, a40826, a41512, with different residues.
* x24601 is singular: it is the only bit that pairs to 4 with anything.

So the checkpoint's message is not unique — there are two siblings with different loaded
constants and identical structure, which is exactly what an invariant-hitting search needs.

## 5. The third channel, tested

The never-explored channel `x34606 = U*(1-V)` is reached by `{x24601}` alone.  Its response
matrix (`s11/data/resp_C.pkl`) gives a 2,027 x 1,470 system, rank 1,470, **also inconsistent** —
but with far smaller certificates than channel A's:

    {a688, a40608}                              (two rows, both failing)
    {a1618, a2423, a31670}
    {a688, a26731, a31672, a33929}
    {a25676, a33792, a40562, a40623, a42245}    (a25676 and a42245 in ONE equation each)
    ... 8 in total

Several members sit in exactly one equation, so a cheap hitting set looked possible — but every
certificate needs a688 (11 equations) or a1618 (14), and the drop search confirms nothing beats
7 there either (`s11/chanC2.py`).

## 6. Where this leaves the search

The obstruction is now an algebraic target rather than a plateau:

    find a message with  INV_0 = ... = INV_5 = 0  in GF(p)

with a known factorisation (inv5 depends on 37 bits; inv0/1/4 on 78 in the V=0 regime), a known
degeneracy (inv5's image is tiny), a fast exact oracle (0.08 s), and three known sibling messages
at the 4-failure minimum.  That is a genuinely different problem from the one every earlier
session was attacking, and it is the one worth attacking next.

## Files
    s11/bits1..4.py    census of the 256 real bits (trees 88/90/37/41, 2 private pins each)
    s11/bits5,7,8.py   exact weight-1 and exhaustive weight-2 message scans
    s11/bits10..16.py  invariants: conservation, dependency, degeneracy, C-subset search
    s11/chanC,chanC2.py the third channel's system and its drop search
