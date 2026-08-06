# Session 11, Part X — the checkpoint's obstruction, solved exactly

Everything below is exact integer algebra on the 39,026 checkpoint
(`best/new_instance_partial_39026.json`), not a search heuristic.  It supersedes the
"rigidity" language of Parts VI–VIII, which measured restricted move sets rather than the
instance (see the Part IX correction).

## 1. The checkpoint has exactly seven nonzero atoms, and they span exactly seven equations

    broken gate atoms   22229 22230 35758 35761 35762
    broken check atoms  35759 35760
    failing equations   12231 12270 12350 14584 18673 22044 29125

Every other one of the 42,267 atoms is exactly zero, so the whole defect lives in nine atoms
(the seven above plus 35756, 35757, which are zero but share the same equations).  Fifteen
equations touch those nine atoms; eight hold, seven fail.

## 2. The nine atoms, decoded

    z0 = a35756 = x1844   - p*x21574         z1 = a35757 = S - x1844
    y3 = a35758 = x29854  - p*x1329          y4 = a35759 = T - x29854
    y5 = a35760 = x31864  - p*x10903         y6 = a35761 = U + x31864
    y7 = a35762 = x642    - p*x17325         y1 = a22229 = A - 7376877*x642
    y2 = a22230 = x28730  - p*x9413

    S = x1956*x17065     T = 5113045*x7075*x9118     U = x7075*x8731
    A = x7068 - x2099    B = x28730

Nine variables — x1844 x21574 x642 x9413 x1329 x29854 x10903 x31864 x17325 — occur in **no
atom outside these fifteen equations**.  They are private knobs.

## 3. The reachable set is an explicit coset (`s11/local1.py`)

Eliminating the nine knobs, the attainable (z0,z1,y3..y7,y1,y2) are exactly

    z0 + z1        == S   (mod p)
    y3 + y4        == T   (mod p)
    y6 - y5        == U   (mod p)
    y2             == B   (mod p)
    y1 + 7376877*y7 == A  (mod 7376877*p)

and nothing further.  Each of the fifteen equations is a linear form in the nine values with
zero constant term.  So "how many equations must fail here" is a finite integer feasibility
question.  Enumerating drop-sets by size:

    drop 0..6 : infeasible          drop 7 : feasible

**Seven is exactly optimal in this neighbourhood** — proved, not observed.

## 4. What the rest of the circuit has to deliver

The 15x9 coefficient matrix has rank 9, so all fifteen equations hold only when all nine atoms
vanish, which by §3 needs

    p | x9118        p | x8731        p | x28730        7376877*p | (x7068 - x2099)

The fourth is not independent: x7068 = x2099 + 7376877*x642 is a gate and x642 = p*x17325 is a
gate, so it follows from the other three plus the gates.  **Three congruences are the entire
obstruction at the checkpoint.**

## 5. They are reachable — x9118 and x8731 are FREE variables

This is the point every earlier session missed.  x9118 and x8731 are inputs, not computed
values; x28730 is defined by the very gate (a22230) that demands its divisibility.  So round
x9118 and x8731 down to multiples of p, set the four handles to what the checks want, and
ripple the gates (`s11/fix2.py`):

    checkpoint      7 broken atoms  ->  7 failing equations   (score 39,026)
    after fix2      4 broken atoms  -> 29 failing equations   (score 39,004)

All nine atoms of §2 are now zero.  The checkpoint's obstruction is **cleared**, for the first
time in eleven sessions.  The cost is relocation, not removal.

## 6. What is left, and why

The four survivors are all checks of one shape — *a free input variable must equal a computed
copy*:

    a7930  = 9367949*(x24548 - x25442) + ...      x24548 FREE
    a29539 = 12846437*(x14853 - x1308) - x29967   x14853 FREE
    a40826 = 25692874*(x14853 - x1308)            (1 equation)
    a41512 = -252934623*(x24548 - x25442) + ...   (1 equation)

x642 had to move to a multiple of p, so x7068 = x2099 + 7376877*x642 moved, and x7068 is
copied all over the circuit.  Every mirror is free and *can* follow, but each mirror also feeds
gates that feed further checks, so the repair fans out.  Setting the mirrors one at a time
(`s11/repair.py`, `repair2.py`, `fixpoint.py`) cascades and stalls at 11 broken atoms.

The honest linear version (`s11/resp.py`, `resp2.py`) measures the exact ripple response of
every candidate free variable — only CHECK atoms need rows, because a ripple keeps every gate
atom at zero automatically.  Two hops out of the fix2 state:

    386 check rows x 401 exactly-affine columns   (221 further columns are booleans, non-affine)
    rationally consistent  ...  but NOT solvable mod p, and not over Z

So the copy network is *rationally* repairable and the obstruction is integrality — the
p-quantised handles again, now pinned to a concrete 386-row system instead of a vague claim.

## 7. Where this leaves the search

Two concrete, score-relevant targets, both now well posed:

1. **Cheapest drop set.**  86 of the 386 rows lie in exactly one equation (390 of 2,417 at
   three hops from the checkpoint).  If the p-obstruction can be pushed onto at most six
   single-equation aggregate checks, the score beats 39,026.  `s11/drop2.py` searches this by
   cost rather than cardinality.
2. **The 221 boolean columns.**  They are the message bits and the only non-affine freedom in
   the system; they are exactly what a rational-but-not-integral system needs.

## Files
    s11/diag1..5.py   the checkpoint's defect, decoded
    s11/local1.py     the exact coset + drop-set enumeration  (proves 7 optimal locally)
    s11/local2.py     rank 9 => the three congruences
    s11/gp1..3.py     global mod-p census; p-handle census
    s11/fix1,fix2.py  clearing the obstruction (4 broken atoms)
    s11/fix7.py       mirror repair costs
    s11/repair,repair2,fixpoint.py   cascade experiments
    s11/resp,resp2.py exact ripple-response system + diagnosis
    s11/drop2.py      cheapest-drop-set search
