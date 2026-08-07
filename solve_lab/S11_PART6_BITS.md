# Session 11, Part XVII — the obstruction as conserved quantities of the message

Part XIV reduced the instance to a GF(p) circuit and showed no continuous move clears it.  This
part turns that negative into a usable positive: the obstruction certificates are **invariants**,
computable from any message in 0.08 s, and they factor over the OR-trees.

## 1. Certificates are conserved quantities

A certificate y satisfies `y . J = 0` over every continuous knob, so

        INV_y  =  sum_a  y_a * r_a          (r = check residues)

is **constant under every continuous move**.  At the state where y was derived this is true by
construction; the question that matters is whether it survives a change of message.

Two tests, and the second is the honest one.  `s11/bits10.py` found 39 of 40 random perturbations
leaving the invariants unchanged — but that is mostly vacuous, because most knobs do not touch
the certificate rows at all.  `s11/bits19.py` redoes it properly, sampling 150 *live* knobs at the
exact derivation state and at sibling messages reached by swapping the C bit and its two loads:

    gmp16_base (where y was derived) : 2 knobs have an effect; all 6 certificates annihilate both
    sibling x2081 -> x4287           : 1 knob  has an effect; all 6 annihilate it
    sibling x2081 -> x13195          : 1 knob  has an effect; all 6 annihilate it

So the certificates do carry to neighbouring messages.  The evidence is thin (only 1-2 live knobs
touch those rows at all), so treat cross-message validity as *supported, not proved*.

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

Exhaustive over C-subsets of weight <= 3 (`s11/bits16.py`) and random over all weights
(`s11/bits17.py`):

    8,474 subsets of weight <= 3 :  338 distinct values;  one value hit 2,672 times
    1,600 random, all weights    :   82 distinct values;  one value hit 1,063 times
    per weight class             :  ~15-28 distinct values each
    zero: never, in ~10,000 evaluations

So inv5 is roughly 25-to-1 degenerate — nothing like a random GF(p) function.  That is a real
structural handle, and it also means the sampled image is small enough that the absence of zero is
informative rather than a needle-in-a-haystack artefact — though 10,000 of 2^37 subsets is still
a small sample, so it is a signal, not a proof.

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

---

# Part XVIII — invariant 5, enumerated completely

## 7. It sees only 18 of the 256 bits

Classifying the 37 C bits by the value of inv5 on the weight-1 message (`s11/bits20.py`):

    19 bits give the SAME value as the empty C-subset  -- inert for inv5
    18 bits each give a distinct value                 -- active

So inv5 is a function of a subset of **18** bits.  It is *not* additive on them (0 of 40 random
subsets matched the sum of singles, `s11/bits21.py`) — consistent with OR-tree saturation, where
what matters is which subtree fires rather than how many bits fire.

## 8. A cone evaluator, and the complete enumeration

inv5 needs only the transitive input cone of its 8 atoms: **2,888 variables of 38,748**, 2,455 SCC
components of 30,575.  Evaluating just that cone takes **3.9 ms** instead of 80 ms, and it matches
the full evaluator exactly on every test (`s11/bits23.py`).  That makes 2^18 tractable:

    s11/bits24.py:  262,144 C-subsets enumerated EXHAUSTIVELY  (19 min)
                    232 distinct values of inv5
                    zeros: 0

262,144 inputs collapsing onto 232 outputs is extreme degeneracy — multiplicity profile
`[(1,83),(3,54),(127,21),(381,18),(57,11),(174,9)]`, i.e. 83 values hit exactly once and a handful
hit hundreds of times.  **Zero is not in the image.**

## 9. What that does and does not prove

The scoping matters, and a full knob sweep settles it (`s11/bits25.py`, all 1,470 live non-bit
knobs, not the thin 150-sample of bits19):

    checkpoint {x2081,x24601} : 14 knobs touch the certificate rows; all 6 certificates
                                annihilate all 14   -> INV is EXACTLY conserved here
    sibling  {x4287,x24601}   : 14 knobs touch;      annihilate 12 of 14 (13 for inv5)
                                -> NOT exactly conserved there

So:

* **At the checkpoint's own message the result is a proof** (within the affine model): INV_5 is
  exactly conserved, it is nonzero, and it must vanish for a full mod-p solution.  That message
  cannot be completed — which is finally a *reason* for the eleven-session plateau at 39,026,
  rather than another failed search.
* **Across the other 2^18 messages it is a screen, not a proof.** The true certificate drifts with
  the message (12/14 rather than 14/14), so the enumerated quantity is the checkpoint's functional
  evaluated elsewhere. Never hitting zero in 262,144 tries is strong evidence that channel A is
  dead, and no more than that.

The honest reading: channel A is very likely unable to produce a full solution for *any* message,
and the effort belongs in channels B (`x5647`, 3 failing checks) and C (`x34606`, certificates as
small as two rows).  Deriving the exact certificate for a shortlist of messages there — 15 min
each with the cached pipeline — is the concrete next move.

---

# Part XIX — the configuration space, swept completely

## 10. Fifteen classes, and why the checkpoint's message is special

`U = OR(A-leaf, B-leaf)`, `V = OR(C-leaf, D-leaf)`, and the two mirrors are the ANDs:
`x38170 = A-leaf AND B-leaf`, `x3896 = C-leaf AND D-leaf`.  So the configuration space is
(which of A,B fire) x (which of C,D fire) = 16 classes, one killed by the OR gate.  `s11/cls1.py`
builds a representative message for each; `s11/cls2.py` repeats it with the *best* bit from each
tree (from the weight-1 scan):

    U-side     V-side     m2 m1  best message              failing checks
    ('B',)     ('C',)     0  0   {x2081, x24601}              4
    ('A','B')  ('C',)     1  0   {x1413, x2081, x24601}      14
    ('B',)     ('C','D')  0  1   {x2081, x3545, x24601}      14
    ('B',)     ()         0  0   {x24601}                    15
    ...
    ('A','B')  ('C','D')  1  1   {x1413,x2081,x3545,x24601}  35

Two clean regularities: **the count grows with the number of trees fired**, and **turning a
mirror on always costs** — which makes sense, since a mirror being 1 activates the checks it
gates (a26719/a26721/a26723 for x3896; the a688/a1618/a40608 core for x38170).  With the mirror
off those checks read `0 = p*handle` and are vacuous.

The checkpoint's class — one B bit, one C bit, both mirrors off — is uniquely good, and within it
only three bit pairs reach 4.  That is why eleven sessions converged there.

## 11. Solving a channel, not just scoring it

Raw counts are biased: every message inherits whatever free-input background it was built on, and
`gmp16_base` was tuned for channel A.  The fair measure is the residual *after* solving the
continuous system (`s11/solveres.py`):

    channel C ({x24601}, U=1,V=0):  15 failing checks  ->  3 after one linear solve
                                    residual a2423, a10506, a26731  (27 equations)

For channel A the same one-shot solve makes things worse (6 -> 15), which is a nonlinearity
artefact — the required delta is large and the affine model stops holding, so that number is not
channel A's residual and should not be read as one.  The exact statements per channel are the
certificates, not the applied solves.

Channel C's residual atoms (a2423, a10506, a26731) are **different** from channel A's
(a7930, a29539, a40826, a41512).  The obstruction genuinely changes shape between channels, which
is precisely why channel A's invariant enumeration does not settle the instance.

## 12. All three channels, side by side

Each channel's response matrix was built from its own properly forward-evaluated state
(`s11/data/resp_modp.pkl`, `resp_B.pkl`, `resp_C.pkl`), and each has its own certificates:

    channel      message        failing checks   obstruction directions   residual atoms
    A  x15298    {2081,24601}         4                   6               a7930 a29539 a40826 a41512
    B  x5647     {91,490}             3                   4               a26719 a26721 (a26723 is NOT obstructed)
    C  x34606    {24601}             15                   8               a2423 a10506 a26731  (after one solve)

**Channel B is the closest to solvable**: fewest failing checks, fewest obstruction directions,
and one of its three failures — a26723 — appears in no certificate at all, so it is fixable.
Its four certificates all pin a26719 and a26721 together, tied to a688, a1618, a21050, a29253,
a29539, a38567 and a40065.

Its invariants do not factor as sharply as channel A's inv5.  Within channel B (U = 0, so no A or
B bits) they are moved by **42 V-side bits** (22 C, 20 D) — 2^42, beyond exhaustive reach — and
they show the same signature as channel A's (`s11/chanB3.py`, 700 channel-B messages):

    invB0: 128 distinct values, zero never;  invB1: 325, never
    invB2: 143 distinct values, zero never;  invB3: 254, never
    number of invariants simultaneously zero: 0 in all 700

## 13. The consistent signature, and the honest bottom line

Every channel tells the same story: a handful of conserved quantities, each far more degenerate
than a random GF(p) function, none of which is ever observed to vanish — exhaustively over 2^18
in channel A, by sampling in B and C.

What is proved: **at the checkpoint's own message, the obstruction is exact** — the certificates
annihilate every knob that touches them, so that message cannot be completed to a full mod-p
solution, and 39,026 is the ceiling there.  That is the first genuine explanation of the plateau.

What is not proved: that no message works.  The certificates drift between messages (12 of 14
knobs annihilated at a sibling), so the cross-message enumerations are strong screens rather than
proofs, and the instance was presumably generated from a witness that satisfies everything.

The remaining move, precisely stated: pick a shortlist of channel-B messages, derive the exact
certificate for each (13 min per message with the cached pipeline), and test **that** message's
own invariants for vanishing — rather than reusing another message's functional.  Channel B is the
right place because it has the fewest obstruction directions and one already-free failure.
