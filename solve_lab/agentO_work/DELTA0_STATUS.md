# READ THIS BEFORE USING `DELTA0_FOR_M.json` — the target is REACHABLE-IN-LATTICE but BLOCKED

Agent O.  I emitted `DELTA0_FOR_M.json` / `.md` before I had priced the two open carriers in
frame B.  I have now priced them.  **δ₀ is a valid lattice target and it is not realisable over
the knob set that carries it.**  The handoff is still worth having — the shift vector and the
region description are correct and reusable — but do not spend time trying to realise δ₀ as
written without reading this first.

## What I measured (frame B, `frameB.Frame([642,28730,29854,31864])`)
Frame B reproduces the witness bit-for-bit: score 39,026, failing
`[12231,12270,12350,14584,18673,22044,29125]`, **0 variables differing**.

Both open carriers turned out to be *free inputs* in frame B, so I redid the region solve
natively there over 12 knobs:
`{642, 1329, 9413, 10903, 17325, 28730, 29854, 31864}` (region-private) `+ {7068, 4432, 8731, 9118}` (carriers).
Those 12 knobs reach exactly **12 check atoms and 29 equations**, and **all 7 witness failures
are inside that set** — nothing is unreachable.

## The result: a 1-for-1 trade, seven ways
Every one of the 7 failing equations can be bought individually — and **every purchase costs
exactly `eq8680`**.  Score stays 39,026 in all seven cases; the failing set merely rotates.
No subset of size ≥ 2 is buyable.

## Why — and it is one atom, not a lattice accident
`eq8680` is the only equation containing atom **`a37887`**, and `a37887` is a **perfect square**:
its source is literally `(S) * (S)`.  At the witness `S = 0`.  So "eq8680 holds" is not a
quadratic obstruction — it is the **linear** constraint `S = 0`, and

    dS/dx_4432 = +1 ,   dS/dx_28730 = -1 ,   dS/d(every other knob) = 0.

So `S = 0` is exactly `δx_4432 = δx_28730`.  But `x_4432` is the *sole* carrier of the `a23618`
("L") shift and `x_28730` is the private handle already in the region.  **`S = 0` collapses the
L direction onto the x_28730 direction, annihilating precisely the new degree of freedom δ₀
needs.**  With `S = 0` added as an explicit row, **nothing is buyable at all** — not one of the
seven.

My earlier atom-level model could not see this for two reasons, both now understood: it drops
`a37887` as nonlinear (it is quadratic), and `a37887` is the one atom *outside* the region that
ties the two carriers together.

## Status of each of the four shifts
| shift | carrier | status |
|---|---|---|
| `a36662` = `x_7075*x_8731` | `x_8731` | free (zero collateral in frame B) |
| `a36660` = `5113045*x_9118` | `x_9118` | free (zero collateral in frame B) |
| `a23616` = `x_7068 - x_2099` | `x_7068` | free input in frame B; only matters mod 7376877 |
| `a23618` = `x_4432 - x_19964` | `x_4432` | **BLOCKED by `S = 0`** — this is the binding one |

## Scoped claim
Over those 12 frame-B knobs, with every other free input at the witness's values, **39,026 is
exactly optimal**, and the binding constraint is `S = 0` (equivalently `eq8680`, equivalently
`a37887`).  I claim nothing outside that knob set.

This independently confirms and sharpens agent H's "a22231 buys 1 row and costs eq8680,
exactly": it holds for **all seven** failing equations, and the mechanism is `S = 0` forcing
`δx_4432 = δx_28730`.

## I then searched for an independent S-mover — exhaustively, and it is negative
`a37887` is supported by **26 free inputs**; **17** of them actually move `S`.  Two are
`x_4432` (dS = +1) and `x_28730` (dS = −1), both already in the 12 and both zero-collateral.
The other **15** were each added to the knob set and the exact maxsat re-run:

| knob | dS | outside-equation collateral | result |
|---|---|---|---|
| x_6947 | 260b | 10 | nothing buyable |
| x_33168 | 260b | 11 | nothing buyable |
| x_950, x_8976, x_35531, x_22526 | — | 13 | nothing buyable |
| x_12553, x_10422, x_34600, x_11099, x_15324 | — | 14 | nothing buyable |
| x_3629, x_15120 | — | 15 | nothing buyable |
| x_4287 | 3042b | 29 | nothing buyable |
| x_2081 (the selector) | 295b | 117 | nothing buyable (131 rows) |

**Not one of the seven failing equations becomes buyable with any of them.**  So the block is
not "we lack a carrier for S" — S has 17 carriers and none of them helps.

## Final scoped claim
**39,026 is exactly optimal** over the 12 frame-B knobs
`{642, 1329, 9413, 10903, 17325, 28730, 29854, 31864, 7068, 4432, 8731, 9118}`,
and over every 13-knob extension of it by a free input that moves `S`, with all other free
inputs at the witness's values.  Outside those knob sets I claim nothing.

δ₀ remains a correct lattice target for the region system in isolation; it is the coupling to
`a37887` — invisible to any model that drops quadratic atoms — that makes it unreachable.
