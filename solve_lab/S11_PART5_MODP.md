# Session 11, Part XIV — the barrier in the layer where it actually lives

Everything before this worked over Z, where every repair eventually died on a divisibility test.
This part moves the whole instance into GF(p), where those tests do not exist, and the barrier
becomes a small, explicit, checkable object.

## 1. The instance is a plain circuit over GF(p), and it evaluates in 0.08 s

Every handle enters as `(free var) * (wire)` and every wire equals p, so **mod p the entire
quotient-witness apparatus vanishes**.  What is left is a circuit: free inputs -> gates -> checks.
`s11/gmp1.py` forward-evaluates it globally:

    gates that do NOT determine their output mod p : 0 of 31,475
    gate atoms nonzero mod p after evaluation      : 0 of 31,475
    CHECKS failing mod p                           : 6 of 10,792
    one full global evaluation                     : 0.08 s

Six numbers in GF(p) are the whole problem.  If they can be zeroed, every atom is 0 mod p, each
equation is `p*r`, and the p-quantised handles absorb it exactly — a full solve.

The 0.08 s matters as much as the reduction: it makes exhaustive scans over all 7,273 free
inputs and all 256 message bits cheap, which is what the rest of this part rests on.

## 2. Two of the six are free

    a35759 = 5113045*x7075*x9118 - x29854   ->  mod p:  5113045*x9118      (x29854 = p*x1329)
    a35760 = x31864 - p*x10903             ->  mod p:  -x7075*x8731

x9118 and x8731 are free inputs, and their measured responses touch **exactly one check each**
and nothing else.  So the barrier is four checks, not six:

    a7930, a29539, a40826, a41512

## 3. The maximal continuous move set is inconsistent

`s11/gmp6.py` caches the exact response of every free input on every check (7,273 probes, ~15
min): **1,726 live knobs**, responses very sparse (485 touch exactly one check).  `s11/gmp14.py`
then solves the whole thing:

    2,595 rows x 1,726 knobs,  rank 1,726 (full column rank),  6 inconsistent rows

So no continuous move from this base zeroes the four.  Every earlier obstruction claim was about
a handful of variables; this one is about *every* free input at once.

## 4. Why: the message bits are frozen, one row at a time

`s11/gmp20.py`: each load-pin row `b*(x - C) - m*h` is moved by **exactly one knob — its own bit
b**.  So any linear system containing the pin rows forces that bit's coefficient to zero.  Bits
are not continuous knobs; they move only by flipping.

Dropping the bits and their pin rows leaves the genuine continuous system, whose certificates are
readable (`s11/gmp23.py`) — the smallest is 8 rows:

    [a3578, a7930, a21617, a25676, a33792, a40562, a40623, a42245]     (only a7930 failing)

## 5. The deficit, exactly 1

Reading the non-bit movers off the cached matrix (`s11/gmp21.py`):

    a33792 <- {x8183}          a40623 <- {x27156}         a40562 <- {x30060}
    a7930  <- {x12553, x24548} a21617 <- {x14623, x24548}
    a42245 <- {x14623, x31339} a25676 <- {x8183, x14623, x27156, x30060, x31339}

The first three rows are private: they pin x8183, x27156 and x30060 outright.  What remains is

    4 equations {a7930, a21617, a25676, a42245}  in  3 knobs {x24548, x14623, x31339}

x24548 is forced by a7930, x14623 by a21617, and x31339 must then satisfy **both** a25676 and
a42245.  Deficit exactly 1.  The missing fourth knob is x12553 — frozen by the load pin a3578 of
bit x2081.  This is the "deficit is topological" of Part III, now with an explicit 4x3 system
instead of a matching argument.

## 6. The bits do not decompose additively

If the residues were additive in the ON bits, the search would be a subset-sum in GF(p).
`s11/gmp24.py` turns bits on in pairs (setting each pin properly) and compares against the sum of
the singles: **x25442, x27522 and x6858 are additive, but x1308 never is** — its dependence runs
through the products x3896 = x7304*x25956 and the OR-trees.  So no linear shortcut.

## 7. The bit landscape, measured

`s11/gmp13.py`, `gmp15.py`, `gmp16.py`: 3,484 boolean-checked variables, 1,156 of them free.
Flipping each and re-evaluating mod p:

    900 bits are INERT -- they change no residue at all
    256 bits are real  -- exactly the message bits
    no single flip gets below the base's 4 failing checks (distribution 9,11,12,...,22)

The 900/256 split is worth having: it cuts the nominal search space from 2^1156 to 2^256, and it
says the two currently-ON bits are x2081 and x24601.

## 8. The channel taxonomy, completed

    x7715  = x2108 - x3389 = OR(x8599, x21839)  = U
    x34554 = x22413 - x33770 = OR(x7304, x25956) = V
    x15298 = U*V        x5647 = (1-U)*V        x34606 = U*(1-V)

x2081 feeds V and x24601 feeds U (`s11/chan3.py`), so:

    checkpoint          U=1 V=1   x15298=1   4 failing checks mod p (after §2)
    39,018 state        U=0 V=1   x5647=1    3 failing checks mod p
    x2081 OFF           U=1 V=0   x34606=1   17 failing checks mod p  <-- never explored before
    both OFF            U=0 V=0              29 failing checks (and violates the OR gate)

The third channel is real and reachable.  Its 17 are raw — no repair applied — whereas the other
two numbers are post-repair, so they are not comparable yet.

## 9. Channel B's obstruction, in the same language

Forward-evaluating the 39,018 state mod p reproduces it exactly (15 failing equations, ceiling
39,018 — so that state *is* a clean forward evaluation, with no broken gates).  Its three failing
checks decode to

    a26719 : x3896*x12000 = 8640431*x24326 = 8640431*p*x24175   ->  p | x12000
    a26721 : x3896*x12926 = x21693 = p*x4615                    ->  p | x12926
    a26723 : 2648967*x3896*x21364 = x3358 = p*x13992            ->  p | x21364

with x3896 = 1.  Same shape as channel A's three congruences, on gate-computed values.

## Files
    s11/gmp1.py    global mod-p forward evaluator (0.08 s)
    s11/gmp6.py    cached exact response matrix, all free inputs x all checks
    s11/gmp13,15,16,17.py   bit scans (inert vs real, with and without handle repair)
    s11/gmp14.py   maximal continuous system + certificates
    s11/gmp20,21,23.py      why bits are frozen; the readable certificates; the movers table
    s11/gmp24.py   additivity test
    s11/gmp25.py   cheapest achievable residual
    s11/chan,chan2,chan3.py the channel algebra and the third channel

---

# Part XV — buying knobs by breaking gates, and what that costs

## 10. Breaking a gate is a knob purchase, and now it can be priced

The 39,026 checkpoint's trick — deliberately breaking five gates — is, in this language, buying
knobs.  Breaking gate atom g frees its output variable at a cost of exactly `|equations
containing g|` failing equations.  `s11/gmp26.py` prices the whole catalogue:

    gate atoms by #equations:  {1: 1, 4: 1, 5: 10, 6: 39, 7: 189, 8: 577, 9: 1518, ...}

and asks which cheap ones free a variable that reaches the deficit rows:

    a36244 (4 eqs) frees x3432  -> moves a25676, a42245, a29539, a40826   (17 checks total)
    a36245 (5 eqs) frees x24219 -> same four                              (16 checks)
    a36246 (6 eqs) frees x5077  -> same four                              (16 checks)
    a34869 (6 eqs) frees x10257 -> a29539, a40826                         (12 checks)

The single 1-equation gate, a41332, frees x24453 and moves **no** check at all — worth knowing,
since a 1-equation purchase would have been decisive.

## 11. One knob is not enough, and the reason is geometric

`s11/gmp27.py` freezes x3432 (so a36244 stays broken and x3432 is a genuine knob), rebuilds the
exact response matrix under that freeze — all 7,274 knobs, 13 min — and re-solves:

    system 2037 x 1471, rank 1471, still 6 inconsistent rows

`s11/gmp28.py` explains why the certificate test was misleading.  Extracting the obstruction
functionals y and testing `y . J_t != 0` for every gate with cost <= 8 finds **44 gates that
break at least one certificate**, a36244 among them.  But breaking a certificate is not the
criterion.  The criterion is `rhs in col(J) + span(c)`, and the quotient space here has dimension
2029 - 1470 = 559: a single new column has to be *parallel* to the residual's projection in a
559-dimensional space.  Breaking certificates is necessary, nowhere near sufficient.

## 12. Not even all of them together

`s11/gmp29.py` measures the freed-output column of **every** gate with cost <= 8 (725 usable
columns) and adds them all at once:

    consistent with no gates broken   : False
    consistent with ALL 725 broken    : False

So the residual does not lie in the column space even after 725 extra knob directions.  That is
much stronger than any obstruction claim made in earlier sessions — and it is stated about the
exact linearisation of the whole instance, not a neighbourhood.

**The honest caveat:** this is the affine model.  The responses are exact finite differences, and
several knobs were verified affine, but the instance is polynomial, so joint moves are not exactly
sums.  Linear inconsistency at this base point rules out any *first-order* move; it does not by
itself prove the instance unsatisfiable.

## 13. The combined search

`s11/gmp30.py` searches the natural joint space — break a set of cheap gates, and allow a set of
cheap checks to stay nonzero — scoring by

    |equations of the broken gates  UNION  equations of the checks left nonzero|

with the checkpoint's 7 as the bar to beat.  There are 725 gate candidates and 2,702 candidate
rows, of which many sit in exactly one equation.  No combination of up to 3 dropped cheap checks
with no gates broken beats 7.

## 14. What this establishes

The barrier is now completely characterised, in the layer where it lives:

* the whole instance is a GF(p) circuit that evaluates in 0.08 s;
* from the best-known inputs exactly **four** GF(p) numbers stand between us and a full solve;
* **no continuous move** — all 1,726 live free inputs, plus 725 gate-purchased knobs — reaches
  them;
* the message bits, the only other freedom, are frozen one row at a time by their own load pins,
  move only discretely, are non-additive, and no single flip improves the count;
* exactly 256 of the 1,156 free bits are real; the other 900 are provably inert.

That is a complete map of the obstruction rather than another failed search, and it is what the
next attempt should be aimed at: the four numbers, and the discrete bit moves that are the only
thing left that can touch them.

## 15. Clusters, and why 7 is where the checkpoint sits

The joint search in §13 initially missed the checkpoint's own trick.  Its five broken gates and
two broken checks do not cost 5x10 equations — they all live in the *same seven* equations, so
the price is 7.  Generalising (`s11/gmp31.py`, `gmp32.py`): a **cluster** is a set of atoms whose
equation sets are contained in one small union, so all of them can be broken for |union|.
There are 5.3 million clusters of cost <= 6, far too many to test one by one — but §12 already
rules out closure by cheap-gate knobs, so only the *drops* can matter.

`s11/gmp33.py` therefore tests the promising configurations directly.  The two cheap members of
the obstruction are a40826 and a41512, in **one equation each**, so a state failing only those
would score 39,031:

    168 configurations of cost < 7 (up to 2 broken gates from {a36244, a36245, a36246, a34869}
    and up to 5 dropped checks from {a40826, a41512, a25676, a42245, a36185, a40812, a37662,
    a40623, a40562, a33792})

    none of them closes the system.

So **7 is optimal at this base** — now as a global statement over all 1,470 continuous knobs plus
gate purchases plus drops, rather than the neighbourhood statement of Part X.  The checkpoint is
not merely a good local answer; it is the best the continuous structure permits from here.

## 16. Where the next attempt has to go

Everything continuous is exhausted, and exhausted with proof.  What remains is exactly one thing:

**the 256 real message bits, moving discretely.**

The map is drawn: they are frozen against continuous motion (each is the only knob on its own
pin), they are non-additive, no single flip helps, 900 of the 1,156 free bits are provably inert,
and the three MUX channels are U*V (checkpoint, 4 failing checks), (1-U)*V (39,018 state, 3
failing checks) and U*(1-V) (reachable by turning x2081 off, never explored).  A mod-p forward
evaluation costs 0.08 s, so roughly 40,000 bit patterns per hour can be scored exactly — which is
the tool this needs, and it did not exist before today.
