# The p-wire route, measured (session "wr")

All scores are exact equation counts from `L.failing_eqs`; the headline one was
re-verified with `python3 checker.py`.  Baseline to beat: **39,026** (7 failing).

## Frames used
* `s10/wr_frame.py` — `F_WIRE` = frame3's six detachments **plus** `x_26064`
  (so `a37694` is a check).  The delivered witness is on-manifold there: 39,026.
* wire-free frame (`wr_kernel2.py`) = frame3 + **all 220 wire members detached**,
  so every wire copy atom becomes a check and a non-uniform wire is realisable.

## 1. Uniform wire value w (the w = 1 idea)  — **39,020, checker-verified**
| w | score | nonzero atoms |
|---|---|---|
| p (control) | 39,026 | the 7 residual |
| 1, −1, 2, −2, 3, 0, −p, 2p, p², p±1 | 39,008 → repaired to **39,020** | 2 |

`s10/wr_engine_w1_x7068_39020.json`:
```
[checker] satisfied 39020/39033  (13 failing)
first 13 failing: [8429,11166,11915,12594,23869,25313,26785,31400,32300,36106,36767,37257,37666]
```
**Only two atoms in the whole 42,267-atom instance are nonzero**: `a37694`
(12 equations) and `a39417` (1 equation, eq 11915).  Everything else — both
residual gadgets, `a7930`, `a29539`, `a40826`, `a41512`, all 1,249 handles — is
exactly zero.  So *the w = 1 claim is true*: the wire value frees the entire
residual.  Its price is 13 against a give-up price of 7.

The 13 is irreducible on the uniform branch:
* all 12 equations of `a37694` carry it with a nonzero coefficient
  (−38, 40, −8, −7, 1, 15, …) and no other atom is nonzero → all 12 fail;
* eq 11915 is `217608357 * (a39417)²` with `a39417` the **only** atom in it, and
  on the diagonal `a39417 = −8·(w − p) ≠ 0`.

## 2. Exact accounting for the placement S = {a37694} (task item 3)
Free-atom relaxation (each atom an independent rational — strictly more freedom
than the circuit has), maximise satisfied over region R(S) with `a37694 ≠ 0`:

| S | region | max satisfied | failing ≥ |
|---|---|---|---|
| {a37694} | 12 | 0 | **12** |
| {a37694, a39417} | 13 | 1 | 12 |
| + copy atoms a37691/2/3 | 14 | 4 | 10 |
| + all five a30970–74 | 19 | 9 | 10 |
| + the boolean cluster (13 atoms) | 22 | 12 | 10 |
| closure round 1 (33 atoms) | 44 | 35 | 9 |
| closure round 2 (69 atoms) | 79 | 71 | 8 |

Calibration: the same relaxation on the *current* 7-atom placement returns
region 12 / satisfied 6 / failing 6, where the true answer is 7 — it understates
by ≥ 1.  Growth-counting bound `|R| − |S| + 1` never drops below 8 either.
Over **Z** (the real constraint) it is worse still: exhaustive search over all
196,624 three-knob systems that admit an *integral* solution at wire base 1
(54,734 of them) found a best of **11** broken rows, versus 13 for uniform.

## 3. Partial / kernel deformation (task item 4)
Linear model: `w_u = p + d_u`, 230 wire-identity atoms, **219 rows, rank 217,
kernel dimension 3** (reproduced).

* 0/1 subsets `d = (1−p)·1_T`, 60 random restarts with full local search,
  T forced to contain a residual multiplier: **uniform T = all 220 is optimal at
  12 broken rows**.  Single members cost 34 (x_22665), 58, 60, 70; all-but-root
  costs 50; root alone costs 62; cheapest member overall is x_3915 at 9.
* Row 37257's wire content is the bare pin alone, so `d_root ≠ 0` always breaks it.
* Kernel deformation (0 identity rows broken) — only realisable with the copy
  gates detached; measured **39,010 → 39,011** after the engine.  The damage is
  exactly the atoms the linear model cannot see:
  * **8 quadratic wire-only checks** `a6017, a15935, a17116, a39967, a41117,
    a41158, a41788, a42158` — 1 equation each, all vanish iff the wire is on the
    diagonal;
  * **3 mixed checks** `a39084, a39417, a41278` — 1 equation each;
  * and the residual does **not** dissolve (the multipliers become 325-digit, so
    the congruences get harder, not easier): 12 more equations.
* Best integral 3-knob deformation at wire base 1 (`x_12752, x_13720, x_18306`,
  80-digit offsets): 11 identity rows, eq 11915 saved — but measured **39,010**,
  because leaving the diagonal switches on the 8 quadratic checks.

## 4. Why it fails, precisely
The wire deformation space splits into exactly two regimes:

* **on the diagonal** (all 220 members equal, any common value ≠ p): the 8
  quadratic wire checks and the 3 mixed checks stay zero, the handles become
  granularity-w and the *entire* residual repairs — but the bare pin `a37694`
  breaks in all 12 of its equations and `a39417 = −8(w−p)` breaks eq 11915.
  **Price 13.**
* **off the diagonal** (including the whole 3-dimensional free kernel): the pin
  can be kept at zero and 0 identity rows need break — but all 8 quadratic wire
  checks plus 3 mixed checks fire (11 equations, 1 apiece) *and* the multipliers
  take 78–325-digit values, so the residual congruences do not dissolve
  (+12 equations).  **Price ≥ 11, measured total 22–23.**

There is no third regime: the diagonal is a 1-dimensional line, the pin's row
(37257) is a nonzero functional on it, and the quadratic checks vanish exactly on
it.  Both regimes cost more than the give-up price of 7.

## 5. Score table (all measured this session)
```
deliverable, p-wire                                  39026   (7 failing)
uniform w != p, before repair                        39008
uniform w = 1, greedy repair                         39011
uniform w = 1, engine, root restored to p            39021   (4-atom placement)
uniform w = 1, engine, full repair                   39020   <-- checker-verified
3-dim kernel deformation                             39010 -> 39011
best integral 3-knob deformation at wire base 1      39010 -> 39011
```
Nothing beat 39,026; `s10/wr_best.json` was not written.
