# Session 11, Part II — the channel taxonomy, and the deficit of 2 made topological

Part I (`S11_SEMANTICS.md`) decoded the circuit and solved the core system with a cubic mod p.
Part II inspects the **39,026 checkpoint in those coordinates**, discovers a much cheaper family
of branches, drives one of them to a verified **39,018**, and finally pins the residual
obstruction to two *structural* control collisions. Everything below is reproducible in `s11/`.

Deliverable unchanged: `best/new_instance_partial_39026.json` (39,026/39,033).
Best reached this session, independently verified with `checker.py`:
`s11/data/finish3_named.json` → **39,018 / 39,033** (15 failing).

---

## 1. Reading the checkpoint in the new coordinates (`s11/wit.py`)

    a = x_8599 = 0    b = x_21839 = 1     c = x_7304 = 1    d = x_25956 = 0
    U = 1  V = 1      x_15298 = U*V = 1
    x_38170 = a*b = 0        x_3896 = c*d = 0        <-- BOTH mirror gates OFF
    message bits on: 2  (x_2081, x_24601)
    x_12186 - x_1308 == 0 (mod p)     x_24908 - x_19083 == 0 (mod p)

> The checkpoint sits in channel **U=V=1 with both mirror gates off**. My Part-I 4-bit
> configuration turned on all of a,b,c,d and lit **both** mirror cores for nothing.
> The checkpoint also fails at a *different* place than my frame: its seven nonzero atoms are
> the `x_2099` arithmetic ladder (22229, 22230, 35758..35762), while `a688/a1618/a40608` are
> already zero there. In my frames it is the other way round.

## 2. The channel taxonomy (`s11/cfg.py`, `s11/two.py`, `s11/uv01.py`)

`a23000 = (1-U)(1-V)` forces U=1 or V=1, and the three data channels are `U*V`,
`(1-U)V`, `U(1-V)`. What actually matters is which of the two mirror gates lights:

| channel | `x_15298` | mirror gates | first core | arithmetic slots |
|---|---|---|---|---|
| U=V=1, `ab=cd=0` (2 bits, e.g. 542+438) | 1 | both off | **active** (2 conditions, 1 knob) | both free (`x_30213`, `x_22162`) |
| U=1,V=0, `ab=1`   (542+47)              | 0 | group-2 on | **dead** | `x_16742` free, `x_12186` via `x_5096` |
| **U=0,V=1, `cd=1`** (490+91)            | 0 | group-1 on | **dead** | `x_24908` via `x_19750`, `x_13682 = x_14853` free |

In every 2-bit configuration **both gaps and all four core quantities are already zero**
(`s11/two.py`) — no cubic needed there at all.

**The best branch is U=0, V=1 with bits (490, 91)**: the group-1 mirror is satisfied from the
start, and after loading the two bits' constants only **8 checks** are nonzero.

## 3. Driving that branch (`s11/uv01build.py`, `s11/tri7.py`)

    x_37892 = x_24908  <- x_19750     drives a688
    x_13682 = x_14853  (free)         drives a1618
    x_1308  <- x_14515                drives a29539
    x_16742 := x_19083                drives a26731

    => a688 = a1618 = a40608 = 0 EXACTLY, mirror zero, both gaps zero, only 5 checks left.

Then the four *linking* checks close against their own free variables:

    a7881 <- x_2751 ,  a21050 <- x_16441 ,  a26839 <- x_18751 ,  a40065 <- x_28955

## 4. The extra divisibility, solved (`s11/quad8640431.py`, `s11/quad3.py`)

The group-1 mirror trio needs more than `x_3719 ≡ x_25118 ≡ 0 (mod p)`:

    a26721 handle x_4615  : delta = -p            -> needs  p | x_12926
    a26723 handle x_13992 : delta = -p            -> needs  p | x_21364
    a26719 handle x_24175 : delta = -8640431*p    -> needs  8640431*p | x_12000

With `x_12000 = 9974121*x_3719 + 15683097*x_25118` this is a condition on `x_12000/p` modulo
**8640431 = 53 x 163027**. Shifting `x_31339 += k*p` and `x_33708 += l*p` preserves the mod-p
mirror, and `gamma(k,l) = x_12000/p` is a bivariate polynomial of **bidegree (2,3)** — verified
by exact interpolation on a grid *and* on a held-out point. Solving it mod 53 and mod 163027
(quadratic in k for each l) and combining by CRT gives **gamma = 0** exactly.

> That closes `a26719/a26721/a26723`. This is a genuine two-stage CRT solve on top of the
> mod-p system; nothing like it was tried in sessions 1-10.

## 5. THE DEFICIT OF 2, LOCALISED — and it is topological

With gamma held at 0 and the mirror trio closed, exactly **two** checks remain
(`s11/closehit2.py`): `a14445` and `a27139`. An exhaustive scan of **all 7,253 non-locked free
inputs** (`s11/last4.py`) against the six live residuals returns:

| residual | # controls | **non-bit** controls |
|---|---|---|
| `a14445`  `x_33129 ≡ x_3757`  | 14 | **[33129]** |
| `a34580`  `x_33708 ≡ x_10170` | 19 | **[33129]** |
| `a27139`  `x_37088 ≡ x_13585` | 5  | **[37088]** |
| `a33796`  `x_31339 ≡ x_6858`  | 32 | **[37088]** |
| mirror `x_3719`, `x_25118`    | 61 | **[]** (message bits only) |
| gamma                          | 66 | [257, 20414, 26489, 29261, 33986] |

> **Two pairs of constraints, each with exactly ONE shared non-bit control.**
> `a14445` and `a34580` are both driven solely by `x_33129`; `a27139` and `a33796` both solely
> by `x_37088`. That is a deficit of exactly **2**.

It is **structural, not numeric**: `x_33129` is the free variable of `a14445` *and* feeds
`x_15111 -> x_20541 -> x_10170`, which is the other side of `a34580`. Changing which message
bits are on changes the *values* but not this topology, so no bit choice removes the collision.

    x_10170 := x_20541 + x_29617        x_15111 := x_33129 * x_35985

## 6. What the deficit costs, and why 15 is the floor here

The two deficits must be absorbed by leaving checks broken. The options, priced in equations:

    break a14445 (14) + a27139 (14)                        = 28   (s11/closehit2)
    break a34580 (13) + a33796 (12)                        = 25
    break the mirror trio a26719+a26721+a26723 (overlapping) = 15   <-- cheapest
    (breaking the mirror frees x_31339 / x_33708 to serve a34580 / a33796,
     and x_33129 / x_37088 then close a14445 / a27139)

Measured exactly: **15 failing equations, score 39,018**, verified by `checker.py`.

> So the channel's floor is 15, against the checkpoint's 7. The checkpoint wins not because
> its deficit is smaller — it is also 2 — but because in *its* channel the two deficits are
> absorbed by the `x_2099` arithmetic ladder, whose seven atoms occupy only **7** equations
> in total. Cheapness of the absorbing set, not the size of the deficit, is what decides.

## 7. Consolidated

* Everything upstream of the deficit is now **solved constructively**: channel choice, the
  arithmetic cluster, both gaps, both cores, the cubic, and the 8640431 CRT condition.
* The obstruction is exactly two control collisions, and they are topological.
* The instance's difficulty is therefore *not* a wall of algebra — it is a rank-1 shortfall
  in the control map, duplicated, placed so that the cheapest absorbing set costs 7 equations.

### Reproduce

    cd solve_lab/s11
    python3 wit.py            # the checkpoint in structural coordinates
    python3 two.py            # channel taxonomy
    python3 uv01build.py      # U=0,V=1 structural solve -> 5 bad
    python3 quad3.py          # bidegree-(2,3) interpolation + CRT -> gamma = 0
    python3 closehit2.py      # -> 2 bad
    python3 last4.py          # exhaustive control scan -> the two collisions
    cd .. && python3 checker.py s11/data/finish3_named.json   # 39018/39033

---

## 8. Part III — the deficit is exactly 2, PROVED by bipartite matching

Everything above is constructive. This section proves the obstruction rather than observing it.

**Bits are genuinely boolean.** `s11/boolform.py` finds explicit checks
`a1430 = x_490² − x_490`, `a1431 = x_12095² − x_12095`, … — one for **every** one of the 256
message bits, each living in 13–14 equations. So the bits are *not* continuous controls; a
non-boolean bit costs 13–14 equations. (This closes the one loophole that would have removed
the deficit outright.)

**The constraint/control bipartite graph** (`s11/hall.py`), built from the exhaustive scan of
all 7,273 free inputs with the bits removed:

    a688   -> [19750]              a1618  -> [14853]            a29539 -> [14515]
    a26731 -> [16742]              a7881  -> [2751]             a21050 -> [16441]
    a26839 -> [18751]              a40065 -> [28955]
    a14445 -> [18751, 33129]       a34580 -> [33129, 33708]
    a27139 -> [2751, 37088]        a33796 -> [31339, 37088]
    mirror3719  -> [14515, 16441, 28955, 31339, 33708]
    mirror25118 -> [14515, 16441, 19750, 28955, 31339, 33708]

    constraints = 14      MAXIMUM MATCHING = 12      DEFICIT = 2
    unmatched: mirror3719, mirror25118

**Hall violator** — 9 constraints whose entire control neighbourhood is only 8 variables:

    {a14445, a21050, a26839, a29539, a34580, a40065, a688, mirror3719, mirror25118}
      -> {14515, 16441, 18751, 19750, 28955, 31339, 33129, 33708}

> **This is the trapdoor, stated exactly.** Not "the system looks rigid", not "7 is an
> invariant" — a Hall violator of size 9 over 8 controls in the constraint/control incidence
> graph, verified against every free input in the instance.

**Pricing the absorbers** (`s11/pairprice.py`). Two constraints must be left unsatisfied; the
cost is the number of equations their atoms occupy:

    a40065 10 | a688 11 | a29539 12 | a21050 12 | a26839 12 | a33796 12
    a34580 13 | a1618 14 | a14445 14 | a27139 14 | a26731 16 | a7881 18
    cheapest PAIR by union:  a688 + a1618 = 15      mirror trio = 15

so the branch floor is **15 failing equations** — exactly what was achieved and verified
(39,018). Also checked (`s11/compensate.py`): no atom in the instance has an equation-footprint
proportional to any of these absorbers, so single-atom equation-space compensation is
unavailable; and the constrained equation-space solve over the full 173-equation region with
26 exact-linear handles (`s11/realise3.py`) returns **NONE**.

## 9. Why the 39,026 checkpoint still wins

The deficit is 2 in every channel examined. What differs is the **price of the absorbing set**:

| channel | absorber | cost |
|---|---|---|
| U=0,V=1 bits (490,91) | mirror trio, or a688+a1618 | **15** |
| U=V=1 `ab=cd=0` (the checkpoint) | the `x_2099` ladder, 7 atoms | **7** |

The checkpoint is better not because its deficit is smaller but because its two unsatisfiable
constraints happen to live in atoms occupying only seven equations between them. Sessions 9–10
proved 7 optimal for that placement by exhaustive subset enumeration; this session independently
explains *why* a deficit exists at all.

> **The single most useful next question:** is there a channel whose 2-deficit is absorbed by
> constraints occupying fewer than 7 equations? Cost is driven entirely by the absorber's
> equation footprint, and the instance contains many 1-equation checks — but none of them sits
> in the tight set of any channel examined so far.

## 10. Part IV — the equation-space optimum, checked properly

Session 9's key correction was "an equation is zero iff its *linear combination of atoms* is
zero" — atoms need not vanish. Part IV applies that here, correctly.

**Region analysis** (`s11/region.py`). For each state, S = failing equations and the *knobs* are
the atoms whose entire equation footprint lies inside S (moving them changes nothing outside).
Solving over the knobs as free integers reports **all** of S recoverable — but that is too
generous: a knob is only movable through an actual free handle.

**Realisable knobs** (`s11/private.py`). Of the six knobs of the 15-equation region, only three
have a *private* handle (a free variable occurring in exactly one atom):

    a26719 <- x_24175   step  -8640431 * p
    a26721 <- x_4615    step        -1 * p
    a26723 <- x_13992   step        -1 * p

The other three (`a26720`, `a26722`, `a28437`) are gates — `x_24326 = x_24175 * x_35019` etc. —
and carry no independent freedom. So the honest system is **15 equations in 3 lattice knobs**.

**Result** (`s11/eqopt.py`, `s11/eqopt2.py`, run over every saved state): no assignment of the
lattice knobs satisfies even |S| − 8 of the failing equations, at any state. Since the knob
steps are multiples of p, they can only repair equations whose atom-sum is already ≡ 0 mod p,
and the failing ones are not.

**Also ruled out** (`s11/compensate.py`): not one atom in the instance has an equation footprint
proportional to any candidate absorber, so single-atom compensation is unavailable; and the
constrained solve over the whole 173-equation region with 26 exact-linear handles is
inconsistent (`s11/realise3.py`).

> Equation space offers no escape here. The atom-space floor of 15 for this branch is also the
> equation-space floor.

## 11. Part V — breaking GATES (the checkpoint's actual trick), priced

Running my pipeline's forward evaluator on the checkpoint destroys it: `fw.forward` recomputes
every gate output, and the checkpoint's score depends on **five broken GATE atoms**
(22229, 22230, 35758, 35761, 35762) plus two broken checks. Forward-evaluating it gives 37
failing — exactly session 10's "forward-eval frame". So the whole session-11 pipeline, which
keeps every gate satisfied, had structurally excluded the checkpoint's strategy.

**Breaking a gate frees its output variable as a new control**, at the cost of the equations
that atom occupies. Pricing all gate atoms (`s11/breakgate.py`, `s11/cheapgates.py`):

    817 gate atoms live in <= 8 equations; the cheapest is a41332 [1 equation] -> x_24453
    (the constant pin x_24453 = 97553848...891), then a36244 [4] -> x_3432, a36245 [5] -> x_24219

Scanned **at a solved state** (the derivatives vanish at the raw baseline), twelve of these
cheap gates move the mirror residuals — which had **no** non-bit control at all:

    a41332 [1 eq] -> x_24453 : moves x_3719
    a36244 [4 eq] -> x_3432  : moves x_25118 and x_3719

So breaking those two costs **5 equations** and would, if they supplied two *independent*
directions, give a perfect matching and a score of 39,028.

**They do not.** A joint 6x6 Newton over
`{x_25118, x_3719, a34580, a33796, a14445, a27139}` with controls
`{3432, 24453, 33708, 31339, 33129, 37088}` closes three residuals and stalls on the other
three at **all 12 random starts** (`s11/joint6.py`): the Jacobian is singular, because
`x_24453` and `x_3432` reach the mirror through the same channel as `x_31339`/`x_33708`. The
gate-breaking lever buys at most **one** dimension, not two — so the deficit drops to 1 at a
cost of 5, and the cheapest remaining single absorber (a40065, 10 equations) brings the total
back to 15. No gain.

> The deficit of 2 survives the one attack that was structurally unavailable to the rest of
> this session. Gate-breaking is real, cheap (1-4 equations), and correctly the checkpoint's
> mechanism — but in this branch the freed directions are not independent of the ones already
> in hand.
