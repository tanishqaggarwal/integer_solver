# Constraints on the 256 selector bits alone

Scope: everything `EQUATIONS.txt` says about the 256 boolean selectors
(`solve_lab/anneal/chain.json -> chain_bit_vars`, identical to field 1 of
`solve_lab/pinrec.json`) **without** invoking the accumulator/ECDLP condition.

Everything below was recomputed from `EQUATIONS.txt` in this directory; no
earlier artefact was trusted. Parse -> 39,033 equations, 46,298 distinct atoms
(`parse_cache.py`, 21 s, cached in `cache.pkl`).

---

## Bottom line

| quantity | value |
|---|---|
| forced bits (`b = 0` / `b = 1`) | **0** |
| equal pairs `b_a = b_b` | **0** |
| opposite pairs `b_a + b_b = 1` | **0** |
| mutually-exclusive pairs `b_a*b_b = 0` | **0** |
| independent GF(2) / Z linear relations among the bits | **0** |
| higher-degree pure-bit relations beyond booleanity | **0** |
| cardinality constraints | **1** - `OR(b_0 ... b_255) = 1` |
| **free bits remaining** | **256** |
| search space | `2^256 -> 2^256 - 1` |

**There is no extra structure to exploit.** The equation system says exactly two
things about the selectors on their own: each is boolean, and they are not all
zero. The second is already implied by the accumulator condition (`k*G = T` with
`T != O`, so `k != 0`). Basis for the linear constraint space: **empty**.

---

## 1. Method

`parse_cache.py` parses every line with `ast`, strips the outer nonzero scalar
factor / square (`LHS = c*CORE` or `CORE*CORE`, so `LHS = 0 <=> CORE = 0`), splits
`CORE` at top-level `+` into `coef * atom`, and expands each atom into a
canonical integer polynomial (gcd-reduced, sign-normalised). It caches:

* `atoms` - 46,298 canonical atom polynomials,
* `eq_terms` - per equation, the list of `(coef, atom_id)`,
* `eq_poly` - the **fully expanded exact polynomial of each equation**.

All support/degree statements below are computed from `eq_poly` (the true
equation), not from the atom decomposition, so they hold regardless of any
modelling assumption.

---

## 2. Every atom that mentions a selector - exhaustive classification

`bit_atoms.py`. **1,836** of the 46,298 atoms mention at least one selector.
Every one of them falls into exactly five families:

| count | shape | meaning |
|---:|---|---|
| 256 | `b - b^2` | booleanity, one per selector |
| 256 | `b - w` | copy: names an alias wire |
| 256 | `1 - b - w` | NOT: names the complement wire |
| 512 | `b*w - C*b - s*v` | pin gate `b*(w - C) = s*v`, **exactly 2 per selector** |
| 68 | `sum c_i (x_i - x_i^2)` | sign-mixed sums of booleanity terms |
| 488 | large random Z-combinations (6-45 variables, deg <= 2) | the "combo" atoms |

The 512 pin gates are the *only* channel through which a selector reaches the
arithmetic part of the circuit. Checked against `pinrec.json`: 512 distinct
pins, exactly 2 per bit, **0 collisions** - no two selectors gate the same
target variable with different constants, so no bit is forced and no pair is
excluded that way.

---

## 3. Atoms and equations supported *only* on the selectors

`scan.py`, `verify.py`.

* **259 atoms** have support inside the 256 bits: the 256 booleanity atoms plus
  exactly three multi-bit atoms (write `t(x) = x - x^2`)
  * `atom#41229` (3 bits, eq **16174**) `2t(x_4381) + 25t(x_26078) - 74t(x_34510)`
  * `atom#45129` (4 bits, eq **33304**) `24t(x_6825) + t(x_15064) - 20t(x_18732) + 25t(x_19435)`
  * `atom#42382` (5 bits, eq **20780**) `33t(x_2143) + 36t(x_2527) + 5t(x_4279) + t(x_5148) - 35t(x_5630)`

* **13 equations** have support inside the 256 bits:

  `7153, 11456, 12810, 13027, 13494, 13905, 16174, 16757, 18154, 19501, 20780,
  30248, 33304`

  Every single one is a pure Z-linear combination of booleanity terms - e.g.

  ```
  eq#11456 :  72*t(x_853) - t(x_2229) + 32*t(x_8319) = 0
  eq#13027 :  40*t(x_17711) - t(x_21392) + 9*t(x_33242) = 0
  eq#19501 : -t(x_47) - 18t(x_1438) - 8t(x_1544) + 27t(x_14644)
             - 28t(x_25104) - 30t(x_29327) = 0
  ```

**Verification** (`verify.py`):
* reducing all 259 atoms and all 13 equations modulo the booleanity ideal
  `<b^2 - b>` (multilinearisation) leaves remainder **0** in every case ->
  each carries *no* information beyond booleanity;
* the degree-1 part after multilinearisation is empty for all 13 ->
  **zero linear relations among the bits**;
* independent numeric check: 200 random 0/1 assignments x 48 equations =
  9,600 exact evaluations, **0 nonzero**.

---

## 4. One-step consequences (bits + exactly one other variable)

`onestep.py`, `verify.py`.

* **518 atoms** with support = bits + one other variable. They are
  * 256 copies `b - w`,
  * 256 NOTs `1 - b - w`,
  * 6 booleanity sums.

  The first 512 are **definitions** of the 512 derived wires, not constraints.
  Eliminating them just renames wires; nothing propagates back onto the bits.

* **35 equations** with support = bits + one other variable
  (`2114, 2277, 3672, 4532, 5931, 8350, 12471, 12806, 12845, 12895, 14088,
  15427, 15600, 16771, 17656, 17704, 19127, 21368, 21591, 25462, 25508, 26018,
  26285, 27056, 27144, 28206, 28696, ...`). All 35 are again pure booleanity
  combinations, now over the bits **and** the extra variable - i.e. they only
  say "this other wire is boolean too". Remainder mod the booleanity ideal of
  (bits + {extra}) is **0** for all 35.

Widening further: **364 equations** in total are pure Z-combinations of
booleanity terms, and together they cover exactly **384 variables** - the 256
selectors plus 128 further boolean wires. 384 is precisely the number of leaf
slots of the combine tree (256 selected + 128 hard-zero pads).

---

## 5. The one genuine pure-bit constraint: `OR(all 256 selectors) = 1`

`ortree2.py`. The OR gate is **not** a single atom, which is why earlier passes
missed it. It is a three-atom gadget:

```
sum  :  s - u - v        ->  s = u + v
prod :  u*v - p          ->  p = u*v
diff :  s - p - o        ->  o = s - p = u + v - u*v = u OR v
```

Matching that triple finds **exactly 383 OR gadgets** - a binary tree over 384
leaves (256 selectors + 128 constant-0 pads), matching the 383 combine gadgets
of `CIRCUIT_STRUCTURE.md`. Climbing the tree from the selectors:

* the root wire class covers **all 256 selectors** (leaf-count histogram tops out
  at 256, one class);
* that class is pinned to 1: `x_2300` (alias of the root `x_9274`) by
  **`atom#16567` = `1 - x_2300`**, which occurs in equations
  **1834, 4279, 8426, 11034, 13146, 15860, 22948, 26944, 32260, 38313**.

So the system asserts `OR(b_0 ... b_255) = 1`, i.e. **not all selectors are zero**
(`k != 0`).

**Independent confirmation** (`closure.py`, `trace_zero.py`, `sweep.py`):
propagating from a concrete selector assignment through the primitive atoms,
the *only* bit pattern out of 338 tested (all-zeros, all-ones, all 256 weight-1,
40 weight-2, 40 random) that reaches a contradiction is **all-zeros**, and the
contradiction is exactly the OR root:

```
x_2855 = b_24136 = 0,  x_35303 = b_25292 = 0        (copy atoms 34295, 34297)
x_618  = x_2855 + x_35303                            (atom 31007)
x_12308 = x_2855 * x_35303                           (atom 201)
x_618  = x_12308 + x_38274                           (atom 21570)
   ==> x_38274 = OR(b_24136, b_25292)
... 383 gadgets up the tree ...
root == x_2300 = 1                                   (atom 16567)
```

At all-zeros the tree evaluates the root to 0 while the pin demands 1 - residual
`1` on `atom#201`. Every other tested pattern propagates contradiction-free.

This constraint removes exactly **one** of the `2^256` bit vectors, and it is
already subsumed by the accumulator condition: `T != O` and `T` is not one of the
256 chain points, so `k = 0` was impossible anyway.

---

## 6. How much booleanity is rigorous vs. modelled

`hard_atoms.py`, `booleanity_proof.py`. Worth stating precisely, because the
usual working model ("every atom vanishes") is *sufficient* but not proven
necessary - 39,033 equations are random Z-combinations of 46,298 atoms.

* **6,592 equations** reduce to a single atom, so those 6,592 atoms are *hard*
  (forced by one equation on their own). Among the 1,836 bit-touching atoms,
  **556 are hard** - but they are all "combo"/"booleanity-sum" atoms; **none** of
  the 256 booleanity atoms, 256 copies, 256 NOTs or 512 pin gates is individually
  an entire equation.
* A rigorous booleanity certificate exists anyway for part of the bits: for
  integer `x`, `t(x) = x - x^2 <= 0` with equality iff `x` in `{0,1}`, so a
  *sign-uniform* equation `sum c_i t(x_i) = 0` forces every `x_i` boolean.
  Propagating that (a proven-boolean variable drops out, which can make another
  equation sign-uniform) gives **58 variables proven boolean, 37 of them
  selectors**, from 18 certificate equations, e.g.

  ```
  eq#2114  -> x_1544, x_14644, x_22960, x_25104, x_30673, x_30823
  eq#6894  -> x_2460, x_3545, x_3872, x_6312, x_21112, x_24376, x_28645, x_37147
  eq#16771 -> x_5584, x_16138, x_22140, x_35079
  eq#17656 -> x_10255, x_32872, x_37438
  ```

  The greedy propagation stops at 37/256; a full LP over the 364 booleanity-
  combination equations would decide whether all 384 follow (attempted with
  exact-rational simplex, did not terminate in the time budget). This does not
  change the answer to the question asked - it only affects whether booleanity
  itself needs the atom model.

---

## 7. What is left: the accumulator condition

For completeness, re-checked against `solve_lab/anneal/instance.py`:
the 256 gated constants are 256 **distinct** points forming a verified doubling
chain `P_i = 2^i*G` on `y^2 = x^3 + B` over `p = 2^256 - 2^32 - 977`, and the
target `T` is neither `O` nor any chain point. So the accumulator condition is

```
sum b_i * 2^i * G  =  T      i.e.   k*G = T,  k = sum b_i 2^i in [0, 2^256)
```

a **single** congruence `k = dlog_G(T) (mod n)` with `n` prime. Because
`2^256 - n = 432420386565659656852420866394968145599 ~ 2^128.3`, that congruence
has one solution in range (two only for the ~2^-128 fraction of targets with
`dlog < 2^256 - n`). It fixes the whole 256-bit word at once and does **not**
decompose into per-bit or per-subset conditions - consistent with sections 1-5:
nothing smaller than the full 256-bit search is exposed.

---

## 8. Files

| file | role |
|---|---|
| `parse_cache.py` | pass 1: exact `ast` parse -> `cache.pkl` (atoms, per-equation terms, per-equation expanded polynomial) |
| `scan.py` | pass 2: atoms/equations whose support is inside bits, and bits+1 |
| `onestep.py` | pass 3: the one-step layer in detail |
| `closure.py` | pass 4: forward propagation from a bit vector, exact evaluation of all 39,033 equations |
| `zeros_probe.py`, `trace_zero.py` | passes 5/9: back-trace of the all-zeros contradiction |
| `shapes.py`, `ortree.py`, `boolsim.py` | passes 6-8: gate census, boolean-layer simulation (negative results: no single-atom OR gate, no pinned shallow boolean wire) |
| `bit_atoms.py` | pass 10: exhaustive classification of the 1,836 bit-touching atoms |
| `verify.py` | pass 11: multilinear reduction proof + numeric spot-check |
| `sweep.py` | pass 12: 338-pattern local-consistency sweep over the 40,356 primitive atoms |
| `ortree2.py` | pass 13: the 383 OR gadgets, the tree, the root pin |
| `hard_atoms.py` | pass 14: which atoms are individually forced by one equation |
| `booleanity_proof.py` | pass 15: sign-uniform booleanity certificates |
| `assemble.py` | pass 16: builds `constraints.json` |
| `constraints.json` | the machine-readable inventory (this document's content) |
| `scan_index.json`, `ortree.json`, `booleanity_proof.json`, `sweep.json`, `hard_atoms.json`, `bit_atom_classes.json`, `boolean_layer.json`, `gates.json`, `boolsim.json`, `closure_results.json` | intermediate data |
| `cache.pkl` | 22 MB parse cache (regenerate with `python3 parse_cache.py`) |

Reproduce in order:

```bash
python3 parse_cache.py && python3 scan.py && python3 verify.py \
  && python3 ortree2.py && python3 hard_atoms.py \
  && python3 booleanity_proof.py && python3 sweep.py && python3 assemble.py
```
