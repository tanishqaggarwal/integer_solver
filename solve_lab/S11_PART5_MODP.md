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
