# The last door: budgeted multi-atom compensation inside `S`

Agent O.  `S = 0` is forced (see `EQ8680_LEMMA.md`) and `a23618` enters the affine form `S` at coefficient +1,
so the `L` direction δ₀ needs exists only if some of `S`'s **other terms** move to keep
`S = 0`.  This is the only route by which the 1-for-1 trade could be leveraged.  Here is what
that channel actually costs.

## 1. There is no free compensator
Every one of `S`'s 20 atoms (18 bracketed groups) appears in **10–18 equations** — none is confined to eq8680 the way
`a37887` is in H's bundled parse.  So no compensation is free; all of it is budgeted.

| atom | coef | #eqs | source |
|---|---|---|---|
| a23618 | 1 | 12 | `x_4432 - x_19964 - x_28730` |
| a23619 | 6 | 12 | `x_23754 - x_26874 * x_6947` |
| a23620 | 15 | 11 | `6122989 * (x_21279 * x_2239) - x_23754` |
| a23621 | −21 | 13 | `x_35619 - x_24490 * x_33168` |
| a23622 | −13 | 13 | `x_21279 * x_31731` |
| a23623 | −13 | 13 | `x_35619` |
| a20448 | 25 | 14 | `x_9629 - x_30095 * x_950` |
| a20449 | 1 | 15 | `x_21279 * x_9106 - 13523997 * x_9629` |
| a20450 | 25 | 16 | `x_18253 - x_4339 * x_15120` |
| a20451 | 28 | 17 | `x_13502 * x_3629 - x_18253` |
| a20452 | 1 | 17 | `x_37720 - x_14466 * x_35531` |
| a20453 | −4 | 15 | `9994531 * (x_13502 * x_8976) - x_37720` |
| a11875 | 23 | 16 | `x_23642 - x_8173 * x_10422` |
| a11876 | −5 | 16 | `x_34600 - x_30108` |
| a11877 | −5 | 16 | `x_23642` |
| a11878 | 20 | 17 | `x_37413 - x_34660 * x_11099` |
| a11879 | −27 | 17 | `x_15324 - x_37254 - 8481759 * x_37413` |
| a11880 | 35 | 15 | `x_23822 - x_30754 * x_22526` |
| a11881 | 17 | 17 | `x_16495 * x_6247 - x_23822` |
| a11882 | −14 | 18 | `x_7945 - x_29084 * x_34868` |

Nine of the twenty are *checks* in E's frame (genuinely movable):
`23618, 23620, 23621, 23622, 20449, 20450, 20452, 11875, 11881`; the rest are definitions.

## 2. The channel was already inside the knob set
The equations those atoms disturb are the **region's own** — 2554, 6816, 8124, 9421, 12231,
12270, 12350, 14584, 22044, 29125 — and every free input that can move any component of `S` is
by definition a carrier of `a37887 = S²` (and eq8680 = S⁴).  All 26 carriers were already in the 34-knob set `K`.
So the compensation channel is not a missing knob; it is a **budget** question.

## 3. The budget, stated as what was actually tested
To beat 39,026 with `j` of the seven bought and `b` satisfied rows broken we need `b < j`.

| budget | scope of the test | solves | result |
|---|---|---|---|
| `j=1, b=0` | **complete** (all 7 singles) | 7 | none |
| `j=2, b=0` | **complete** (all 21 pairs) | 21 | none |
| `j=2, b=1` | **complete** — 21 pairs × each of the 168 satisfied rows *and* the `S=0` row | 3,570 (21 s) | **none** |
| `j=3, b≤2` | **14 of 35 triples enumerated completely** (per triple: b=0, all 168 b=1, all C(168,2)=14,028 b=2) | 198,772 (33 min) | **none** |
| `j=4,5, b≤3` | greedy upper bound only | — | greedy drops 25–26 vs needing <4 |

All 21 pairs are individually feasible, so nothing was pruned away at `j=2` — the negative there
is genuine, not vacuous.  Likewise every triple attempted was feasible on its own.

The greedy pass had flagged `[12231, 12270, 12350]` as dropping *exactly* 3 (net zero — and
greedy only upper-bounds the drops, so the true minimum could have been 2, which would have been
39,027).  Enumerated properly: **exhausted at b≤2, none.**

The 14 triples completed:
`[12231,12270,12350] [12231,12270,14584] [12231,12270,18673] [12231,12270,22044]
[12231,12270,29125] [12231,12350,14584] [12231,12350,18673] [12231,12350,22044]
[12231,12350,29125] [12231,14584,18673] [12231,14584,22044] [12231,14584,29125]
[12231,18673,22044] [12231,18673,29125]`
— i.e. **every triple containing eq12231**, plus nothing else.  The 21 triples not containing
eq12231 were not reached.

## 4. Honest statement of scope
This is **not** "exhausted" at `j=3`: 14 of 35 triples were enumerated completely, the other 21
were not reached within the wall-clock cap.  `j≥4` was touched only by greedy.  The scope of the
whole result remains **34 of 8,751 free inputs, frame B's orientation**.

What *is* complete: `j=1, b=0` and `j=2, b≤1`.  So the 1-for-1 trade is proven unleverageable at
budget 2 over `K`, and unleverageable at budget 3 for every triple containing eq12231.
