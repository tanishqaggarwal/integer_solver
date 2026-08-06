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

---

# Part XI — a second, independent route, and the wall both routes hit

## 8. Turning the bit x4287 on annihilates BOTH congruences at once

    a33881 :  x21279 = x9062*x20434 = x4287*x2081        (x2081 = 1 at the checkpoint)
    a36085 :  x7075  = 1 - x21279

so `x4287 = 1  =>  x21279 = 1  =>  x7075 = 0`, and since

    T = 5113045*x7075*x9118        U = x7075*x8731

both hard congruences of §4 simply vanish — no divisibility on x9118 or x8731 at all.  Better
still, `x31033 = x20434*x31822 = 0` then kills `x22542 = x6418*x31033`, and `x10878` is already
zero, so

    x2099 = x37158 + x25297 = x9118*x21279 = x9118          -- a FREE variable

which means `x7068 = x2099 + 7376877*x642` can be steered directly.  Turning the bit on is a
pinned move: it activates the load pins `a3568` (x31861) and `a3570` (x14865).

## 9. The price: the x21279 channel switches three new congruences on

    a22233 : 6122989*x2239*x21279  = x23754 = p*x6947     =>  p | x2239
    a22235 : x21279*x31731         = -x35619 = -p*x33168  =>  p | x31731
    a19088 : x9106*x21279          = 13523997*x9629       =>  13523997*p | x9106

and the circuit makes all three collapse onto two quantities:

    x2239  = 3494591*x27177 + 14240157*x4306
    x31731 = 15964591*x27177 + 13881285*x4306
    x27177 = x17925^2*(x9118 + x31861 + x6418 + x24453) - x27019^2     -- affine in x9118
    x4306  = (x8731 + x14865)*x17925 - x27019*(x31861 - x9118)         -- affine in x8731

Two knobs, two congruences, both affine — so they are **solvable** (`s11/sw6.py` fits each
response at two points and inverts mod p):

    x27177 = 0 (mod p)   x4306 = 0 (mod p)
    x2239  = 0 (mod p)   x31731 = 0 (mod p)   x9106 = 0 (mod p)     -- all achieved

`a22233` and `a22235` are then repaired exactly.  (`a19088` still wants the extra factor
13523997; shifting x9118 and x8731 by multiples of p moves x9106 without disturbing the mod-p
work, so this is an open but small congruence.)

## 10. Both routes end at the same wall

Solving `x27177 = 0 (mod p)` **pins x9118 mod p**, hence pins `x2099 = x9118`, hence pins
`x7068` mod p — to a residue different from the checkpoint's.  So x7068 moves, exactly as it
does in §5 where `p | x642` forces it.  Either way the free "mirror" inputs that hold copies of
x7068 must follow, and repairing them fans out.

Four independent searches — greedy (`repair.py`), lookahead greedy (`repair2.py`), batch
fixed point (`fixpoint.py`), and beam width 24 / depth 40 (`beam.py`) — run from **both**
routes, all terminate at the same three checks:

    a19297  = x11150*x15298 + x4007                  (11 equations)
    a19299  = x15298*x25739 - 6672769*x29804         (13 equations)
    a30984  = 537773*x15298*x37758 - x35605          (14 equations)

**These three contain no free variable at all.**  `x15298 = U*V` is the live MUX channel, and
every variable in them is gate-computed.  There is nothing to absorb the residue.  That — not
"a factor of p" in the abstract — is the actual wall, and it is now a checkable statement about
three named atoms.

Best reached along these routes: 39,013 (fix2 route) and 38,999 (x4287 route); the deliverable
stays at **39,026**.

## 11. Exact dual certificates of the copy network (`s11/cert.py`)

For the ripple-response system at the fix2 state, the mod-p left kernel gives explicit
obstruction certificates.  Any repair must leave at least one row of each broken:

    support {29539}                                  12 equations
    support {7930, 21617}                            25 equations
    support {31938, 31940, 40826}                    15 equations
    support {7938, 7939, 18691, 18694, 21617, 41512} 35 equations

The singleton is an artefact of the linear model — x14853's response is genuinely quadratic, so
its column is dropped.  Re-basing (setting the mirror exactly, then rebuilding) removes that
certificate and produces a new singleton one row further out.  That is the fan-out, measured.

## 12. Leads left open, in order of promise

1. **A different MUX channel.**  The wall is three checks whose only content is `x15298 = U*V`.
   In a channel with `x15298 = 0` they read `x4007 = 0`, `x29804 = 0`, `x35605 = 0` — plausibly
   vacuous.  Session-11 Part II priced the (490,91) branch at 15 under the *old* restricted-move
   framework; it is worth re-pricing with today's exact tools, because the thing that made 7
   unbeatable in this channel (the x7068 copy network) may not exist there.
2. **13523997 | x9106**, using the p-multiples of x9118 and x8731 as free knobs — completes the
   x4287 channel's local repair.
3. **The 221 boolean response columns.**  They are the only non-affine freedom in the copy
   network, and a rationally-consistent-but-not-integral system is exactly what they could fix.

---

# Part XII — the x5647 branch priced with the same exact tools

Lead 1 of §12 checked out structurally (the three wall checks are vacuous there) but does not
pay off numerically.

## 13. `localopt2.py` — the exact local minimum, for any state

Given a state: A = the atoms occurring in the equations that touch the broken atoms;
KNOBS = variables occurring in no atom outside A (so moving them cannot disturb anything else);
OBJ = **every** equation containing any atom of A — not merely the equations touching the broken
atoms, which was the trap: an atom of A can appear in equations far away.  Then enumerate
drop-sets among the currently-failing equations of OBJ.

Validated: on the 39,026 checkpoint it returns **MINIMUM = 7**, with drop set exactly the seven
known failing equations — independently reproducing `local1.py` by a completely different route.

## 14. The 39,018 state (x5647 channel)

    broken atoms 26719 26721 26723 -- all checks, ZERO broken gates
    a26719 = x3896*x12000 - 8640431*x24326      x24326 = x24175*p,  x24175 free
    a26721 = x3896*x12926 - x21693              x21693 = x4615*p,   x4615  free
    a26723 = 2648967*x3896*x21364 - x3358       x3358  = x13992*p,  x13992 free

with `x3896 = 1` (group-1 mirror on), so the branch's whole obstruction is

    8640431*p | x12000        p | x12926        p | x21364

— structurally the same shape as the checkpoint's three congruences, except that x12000, x12926
and x21364 are **gate-computed**, whereas the checkpoint's x9118 and x8731 are free inputs.
That is why this branch is harder, and it shows up in the numbers:

    localopt2 on 39,018 : movable atoms 84 -> objective 53 equations, 26 knobs
                          drop 0,1,2,3,4,5 all infeasible

    cert.py, ripple-response certificates (obstructions any repair must leave broken):
        1 hop  : no columns at all -- the broken atoms contain no free variable
        2 hops : {26719,26721} {26719,26723} {26719,29539} {26719,40826}
        3 hops : 1095 rows x 863 columns, rank 750;
                 {26719,26721} {26719,26723} {1618,26719,29539}

**Every certificate contains a26719**, at every depth.  A minimum hitting set is {a26719}
itself, costing its 11 equations.  So this branch cannot be pushed below about 11 with any of
these move sets — worse than the checkpoint's 7, and consistent with the drop enumeration.

## 15. Where things stand

    checkpoint channel (x15298 = 1) : local optimum PROVED = 7, obstruction = three congruences
                                      on FREE inputs, clearable, but clearing it moves x7068 and
                                      the copy network walls at a19297/a19299/a30984
    x5647 channel      (x15298 = 0) : those wall checks are vacuous, but the branch's own
                                      obstruction sits on GATE-COMPUTED values and prices at ~11

The deliverable stays at **39,026**.  What is new is that both branches are now priced by exact
integer feasibility rather than by search behaviour, and the two obstructions are written down
explicitly — three congruences each, with the decisive difference being whether the divisible
quantity is a free input or a computed one.

## 16. Sharpest open questions

1. Can `x12000`, `x12926`, `x21364` be steered by the free inputs feeding them?  They are
   sums of gate outputs (`x21219 + x22131`, `x31730 + x36524`, `x6711 + x2853`); the 3-hop
   response system already has 863 columns and still cannot, so it needs the *boolean* columns
   or a different bit pattern in that channel.
2. In the checkpoint channel: `13523997 | x9106` completes the x4287 route's local repair
   (shift x9118, x8731 by multiples of p — they move x9106 without touching the mod-p work).
3. The copy network of x7068: it is rationally consistent and fails only integrally.  The
   221 boolean response columns are the only unexplored freedom there.
