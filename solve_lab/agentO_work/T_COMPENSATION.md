# The last door: budgeted multi-atom compensation inside `T`

Agent O.  `T = 0` is forced (see `EQ8680_LEMMA.md`) and `a23618` enters `T` with coefficient +1,
so the `L` direction δ₀ needs exists only if some of `T`'s **other 19 atoms** move to keep
`T = 0`.  This is the only route by which the 1-for-1 trade could be leveraged.  Here is what
that channel actually costs.

## 1. There is no free compensator
Every one of `T`'s 20 atoms appears in **10–18 equations** — none is confined to eq8680 the way
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
12270, 12350, 14584, 22044, 29125 — and every free input that can move any component of `T` is
by definition a carrier of `a37887 = T²`.  All 26 carriers were already in the 34-knob set `K`.
So the compensation channel is not a missing knob; it is a **budget** question.

## 3. The budget, stated as what was actually tested
To beat 39,026 with `j` of the seven bought and `b` satisfied rows broken we need `b < j`.

| budget | scope of the test | result |
|---|---|---|
| `j=1, b=0` | complete (all 7) | none |
| `j=2, b=0` | complete (all 21 pairs) | none |
| `j=2, b=1` | **complete** — 21 pairs × each of the 168 satisfied rows + the `T=0` row, 3,570 solves, 21 s | **none** |
| `j=3, b≤2` | **exhaustive per triple**, 14,198 solves each; run to a wall-clock cap | see below |
| `j=4,5, b≤3` | greedy upper bound only | greedy drops 25–26, far from a gain |

All 21 pairs are individually feasible, so nothing was pruned away at `j=2` — the negative there
is genuine, not vacuous.

The greedy pass had flagged `[12231, 12270, 12350]` as dropping *exactly* 3 (net zero, and greedy
only upper-bounds the drops, so the true minimum could have been 2).  Enumerated properly:
**exhausted at b≤2, none.**

## 4. Honest statement of scope
This is **not** "exhausted" at `j=3`.  It is: every triple listed in `runs/fb_j3.log` as
`b<=2 exhausted` was enumerated completely; the remainder were not reached within the cap.
`j≥4` was only touched by greedy.  The scope of the whole result remains **34 of 8,751 free
inputs, frame B's orientation**.
