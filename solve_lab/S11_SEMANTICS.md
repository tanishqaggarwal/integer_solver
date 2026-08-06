# Session 11 — the circuit's SEMANTICS decoded (supersedes atom-level models where they conflict)

Ten sessions attacked this instance as *atom algebra* (lattices, kernels, certificates,
hitting sets). This session reads it as what it is: **a compiled arithmetic circuit**, and
decodes the actual program. Every claim below is re-derived from `EQUATIONS.txt` this
session and is reproducible from `s11/`.

Deliverable is unchanged: `best/new_instance_partial_39026.json`, re-verified at
**39,026 / 39,033**. Best score reached constructively in the new frame: **39,013**.

---

## 1. A clean frame prior sessions never used

Set **every free input to 0** and forward-evaluate (`s11/fw.py`). Result:

    all free inputs = 0  ->  6 nonzero CHECK atoms, 28 failing equations (39,005)
    bad checks = [688, 1618, 23000, 39067, 40608, 41211]

Every prior session worked in the *witness frame* (the 39,026 assignment's free inputs),
where the residual is 6 nonzero atoms / 37 failing. The zero frame is strictly cleaner and
splits into **two independent clusters**:

* arithmetic — `a688, a1618, a40608`
* boolean    — `a23000, a39067, a41211`

## 2. `a40608` is not an independent constraint

With `W = x_14257 − 8863713·x_18956`:

    a688   = C − W                      (C = a688's constant term)
    a40608 = W² + qX·W + q0   and   qXY = −2·8863713 , qYY = 8863713² , qY = −8863713·qX

`s11/solveW.py` computes the discriminant exactly: it **is** a perfect square and the
quadratic has a **double root equal to exactly the W that `a688` demands**. So

> `a40608 = (W − C)²`. It vanishes iff `a688` does. Earlier price lists counted them as
> two separate obstructions; they are one.

## 3. The boolean cluster is an OR/AND gadget — decoded exactly

    x_7715  = OR(x_8599 , x_21839)  =: U          x_34554 = OR(x_7304 , x_25956) =: V
    x_9274  = 1  identically                      x_23134 = U·V        x_29237 = U + V

    a23000  =  x_9274 + x_23134 − x_29237  =  1 + UV − (U+V)  =  (1−U)(1−V)

> **`a23000` ⟺ U = 1 or V = 1.** At the all-zero frame both are 0, which is the whole
> reason that cluster fails. `a41211` is a perfect square of a linear form; `a39067` is linear.

Each of `x_8599, x_21839, x_7304, x_25956` is an **OR-tree** over 88 / 90 / 37 / 41 free
boolean inputs. Verified exhaustively (`s11/scan.py`): **every single one of those
88/90/37/41 bits, set alone, switches its node to 1.** 256 bits in total — the "message".

## 4. The data path is a three-way MUX

The handle channels are exactly

    x_15298 = U·V        x_5647 = (1−U)·V        x_34606 = U·(1−V)

— precisely one is active — and they route the two arithmetic quantities:

    x_37892 = x_15298·x_30213 + x_5647·x_24908 + x_34606·x_16742
    x_13682 = x_15298·x_22162 + x_5647·x_14853 + x_34606·x_12186

`a688` needs `x_37892 ≡ −W/8863713 (mod p)`; `a1618` needs `x_13682 ≡ −c₀ (mod p)`.

> Channel **U=V=1** is the only one that makes *both* land on free inputs
> (`x_30213`, `x_22162`). Both were then solved **exactly** (`s11/build2.py`):
> `a688 = a1618 = a40608 = 0`.

`x_2099` is a second MUX of the same shape over free inputs `x_9118, x_31861, x_6418`
selected by `x_4287, x_2081`, guarded by *gated load pins* `b·(X − HUGE − c·p·h) = 0`
(`a3568`, `a3576`).

## 5. Both "cores" have the same shape — and each is TWO conditions, not three

**Group 2** (gated by `x_38170 = x_8599·x_21839`), checks `a26733, a28438, a32342`:

    x_21202 = 11598153·x_25614 + 16335423·x_34220
    x_15286 =  3511239·x_25614 +  9767569·x_34220
    x_32453 =  4677103·x_25614 + 15469317·x_34220

Three linear forms in **two** quantities ⇒ rank 2 ⇒ `x_25614 ≡ x_34220 ≡ 0 (mod p)`, with

    x_25614 = x_17576² − x_18123²·x_32629          x_34220 = x_17576·x_33852 − x_16088·x_18123
    x_32629 = −(x_5096 + x_10261 + x_30454 + x_24453)      x_33852 = x_5096 − x_10261
    x_16088 = x_21589 + x_25199                            x_18123 = x_30454 − x_10261

**Group 1** (gated by `x_3896 = x_7304·x_25956`), checks `a26719, a26721, a26723`: identical
structure, rank 2 in `x_3719, x_25118`:

    x_3719  = x_23776·x_3090² − x_4879²            x_25118 = x_3090·x_26196 − x_2401·x_4879
    x_23776 = x_24453 + x_6083 + x_14515 + x_33708         x_3090  = x_6083 − x_33708
    x_2401  = x_33708 − x_14515                            x_26196 = x_19750 + x_20413
    x_4879  = x_31339 − x_20413

This is the exact form session 9 guessed as `u²·(A·c² − B²)`.

## 6. THE CORES ARE NOT QR-OBSTRUCTED — the real obstacle was a CUBIC

Because `x_28746 = x_18123²` is a perfect square and `x_32629` / `x_23776` are **linear in a
free input**, no square root is needed:

* Group 2 — `x_25614 = 0` is **linear in `x_5096`**; then `x_34220 = 0` is **linear in `x_21589`**.
* Group 1 — eliminating `x_4879` between the two equations gives

        x_23776 · x_2401²  =  x_26196²

  which in `y = x_33708 − x_14515` is a **CUBIC**  `y³ + K·y² − q² ≡ 0 (mod p)`.

> **This is why ten sessions of Jacobians, beam searches, Newton steps and null-space
> arguments all reported "rigid".** A cubic root mod p is invisible to every local /
> first-order method; mod p there are no basins for Newton to descend. Solving it needs
> polynomial factorisation over GF(p), which no prior session tried.

`s11/polyroot.py` implements Cantor–Zassenhaus root finding; `s11/solveA.py` interpolates
the cubic through 5 sample points, takes its roots, and back-substitutes. Result **at the
first usable seed, in under a second**:

    ALL SIX STRUCTURAL TARGETS ZERO
    x_3719 ≡ x_25118 ≡ x_25614 ≡ x_34220 ≡ 0 ,  x_12186 ≡ x_1308 ,  x_24908 ≡ x_19083   (mod p)

Both cores, both gaps, simultaneously — the first time this system has been solved.

## 7. The complete control map (exhaustive over all 7,273 free inputs)

`s11/scangen.py` perturbs **every** free input at a generic point:

    x3719   <- 14515, 16441, 22917, 31339, 33708
    x25118  <- 14515, 16441, 19750, 22917, 31339, 33708
    x25614  <- 5096, 13222, 14681, 28486, 38667
    x34220  <- 5096, 13222, 14681, 21589, 28486, 38667
    n-gap   <- 5096, 14515
    m-gap   <- 19750, 21589

6 equations, 12 controls, two decoupled blocks, triangular order
`(x_5096, x_21589) → (x_14515, x_19750) → cubic for x_33708 → x_31339`.
(At the all-zero point the derivatives of `x25118`/`x34220` vanish identically — they are
products — so the scan **must** be done at a generic point. `s11/scanall.py` shows the
base-point scan returning 0 controls for both, which is the trap that stalls naive Newton.)

## 8. WHERE THE TRAPDOOR ACTUALLY LIVES

Solving the structural system is necessary but not sufficient: **every one of those 12
controls also sits in an always-active *linking* check** of the form
`c·(X − Y) = p·handle`, i.e. `X ≡ Y (mod p)`:

    a21050:  x_16441 ≡ x_4920        a34580:  x_33708 ≡ x_10170       a33796:  x_31339 ≡ x_6858

and each partner has **exactly one** live mod-p control (`s11/partners.py`, cone-exhaustive):

    x_4920 <- x_23210        x_10170 <- x_33129        x_6858 <- x_32125

which is in turn pinned — `x_23210` by **`a38567`, a bit-gated load pin**
`x_91·x_23210 = HUGE·x_91 + x_3556`, and `x_33129`, `x_32125` by further linking checks
(`a14445`, `a35374`) that continue the chain.

> **The sharp statement** (refined by §10 — each of these links is individually closable;
> what binds is the *system*). Each core control heads a chain of mod-p identities that
> terminates in a bit-gated load pin. Once all of them must close simultaneously, the mod-p
> value of every core control is driven by the **message bits**. The structural
> system is solvable in the abstract (§6) but its solution is not reachable from any
> message: the two are connected only through a subset-sum over the 256 load constants.
> That — not "7 is an invariant" — is the precise obstruction.

This also explains, from the construction rather than from failed searches, session 10's
observation that randomising free inputs never helps and that the message space collapses.

## 9. Reproduce

    cd solve_lab/s9 && python3 atomize.py && python3 poly.py && python3 gates.py && python3 fwd.py
    cd ../s11
    python3 fw.py -            # clean frame: 6 bad checks, 39,005
    python3 solveW.py          # a40608 = (W - C)^2, double root
    python3 chan2.py           # the three MUX channels
    python3 build2.py          # a688 = a1618 = a40608 = 0 exactly
    python3 scangen.py         # complete control map (~6 min)
    python3 solveA.py          # cubic solve -> ALL SIX TARGETS ZERO
    python3 partners.py        # the pin chain

---

## 10. Addendum — what actually binds, measured

After the structural solve of §6 lands (all six targets exactly zero), the residue is
**15 bad checks**, and the binding constraint is **handle contention**, not algebra:

* `s11/freedom.py` perturbs each of the 12 controls by 1 and asks, for every check that
  breaks, whether some *other* handle can re-close it. The **linking pins are all
  individually closable** — `a21050`, `a34580`, `a33796`, `a26731`, `a29539`, `a15030`,
  `a9193`, `a31938`, `a31940`. So the controls are *not* individually pinned; §8's chain
  bounds their residues only once the whole system must close at once.
* `s11/ordered.py` closes the residue in hard stages (mirror cores → linking pins → load
  pins → 1-equation checks). Each stage closes 7-11 of its atoms, but the handles overlap:
  closing one class consumes the handle another class needs, and free sweeps then diverge
  (15 → 27 bad).
* `close2.py` (strict monotone acceptance) halts immediately at the first step: **no single
  repair reduces the bad count**, which is the signature of a system that must be closed
  simultaneously rather than sequentially.
* `simul.py` (one exact integer solve over all bad checks against every cone handle,
  74-130 controls) reports the joint system **inconsistent**.

> So the deficit is real and it is a *rank* deficit in the handle map, not a missing idea:
> six structural conditions are satisfiable, but the handles that would realise them are
> shared with the load-pin and linking-check system, and the combined system is
> over-determined. This is the same "deficit of 2" that sessions 9 and 10 kept re-measuring
> in different coordinates — now with its mechanism identified.

**Scores in the clean frame this session:** simplest construction 39,013 (`s11/build3.py`);
structural construction 38,927-39,005. The 39,026 checkpoint's defect placement remains
much cheaper in *equation* terms (7 failing vs 20-28), so it stays the deliverable.
