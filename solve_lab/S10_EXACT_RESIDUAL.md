# Session 10 — the residual in exact closed form, and a proof that 39,026 is optimal here

Deliverable unchanged: **39,026 / 39,033** (`best/new_instance_partial_39026.json`),
re-verified this session with `checker.py` (`satisfied 39026/39033`, failing
`[12231, 12270, 12350, 14584, 18673, 22044, 29125]`) and with a freshly rebuilt
atom model (`s9/atomize.py`: 0 mismatches over all 39,033 raw equations).

Everything below was re-derived from the file, not inherited from earlier sessions.
Tools: `s10/`.

---

## 0. A correction to Session 9's notes

`S9_STRUCTURE.md` §15.1 states that at the 39,026 witness only two atoms are
nonzero and that `D2 ≡ 0 (mod p)` already. That is not the case for the delivered
39,026 assignment. Measured (`s10/resid.py`):

> **7 nonzero atoms: 22229, 22230, 35758, 35759, 35760, 35761, 35762** — living in
> **exactly 12 equations**, of which 5 are satisfied and 7 fail.

The 3-atom state described there is a *different* (39,024-era) point. The rest of
this document uses the measured 7-atom picture.

---

## 1. The residual sub-system

The 7 nonzero atoms occur in exactly these 12 equations:

```
2554  6816  8124  9123  9421 | 12231 12270 12350 14584 18673 22044 29125
     (satisfied)             |              (failing)
```

Every other one of the 39,033 equations has all its atoms zero, so score
= `39033 − 12 + |S|` where `S` is the satisfied subset. The 12×7 coefficient
matrix has **rank 7 over ℚ** (`s10/subsystem.py`), so with no other atom moving,
all 12 hold **iff all 7 atoms vanish**.

Atom sources (all seven):

```
a22229 = x_7068  − x_2099 − 7376877·x_642
a22230 = x_28730 − p·x_9413                     (x_17499 = p)
a35758 = x_29854 − p·x_1329                     (x_22665 = p)
a35759 = 5113045·x_7075·x_9118 − x_29854
a35760 = x_31864 − p·x_10903                    (x_28961 = p)
a35761 = x_7075·x_8731 + x_31864
a35762 = x_642   − p·x_17325                    (x_28599 = p)
```

with `p = 2^256 − 2^32 − 977` and `x_7075 = 1 − x_21279 = 1` on this branch.

---

## 2. The achievable set, in closed form

Writing `A = (a22229, a22230, a35758, a35759, a35760, a35761, a35762)`,
`D = x_7068 − x_2099`, `K2 = x_28730 mod p`, the knobs realise **exactly**

> **(1)  `A1 + 7376877·A7 ≡ D  (mod 7376877·p)`**
> **(2)  `A2 ≡ K2  (mod p)`**
> **(3)  `A3, A4, A5, A6` free**

*Derivation.* `A1 = D − 7376877·t₁` and `A7 = t₁ − p·t₂` give
`A1 + 7376877·A7 = D − 7376877·p·t₂`, which is (1); conversely (1) forces
`7376877 | (D−A1)` and `p | (t₁−A7)`, so both handles are integral.
`A3 = t₄ − p·t₅`, `A4 = 5113045·t₉ − t₄` with `gcd(5113045, p) = 1` make the pair
`(A3, A4)` unconstrained; likewise `(A5, A6)` from `t₆, t₇, t₁₀`.
`A2 = t₈ − p·t₃` would be unconstrained too, but `t₈ = x_28730` is **not** free.

Both congruences were checked numerically at the witness — residues 0 and 0
(`s10/lattice2.py`).

### Which knobs are genuinely free (`s10/isolate.py`)

| knob | ripple | collateral atoms |
|---|---|---|
| `x_642, x_17325, x_9413, x_1329, x_10903, x_29854, x_31864` | 1 var | none |
| `x_9118` | 18 vars | none |
| `x_8731` | 10 vars | none |
| **`x_28730`** | 17 vars | **breaks 7930 and 41512** |

`x_28730` drags `x_4432` (to hold `a22231 = 0`), and `x_4432` feeds atom 7930
(15 equations) and 41512. That is what makes congruence (2) binding.

---

## 3. New: `D` is free modulo p — the 7376877 factor is illusory

`x_7068` can be shifted by any multiple of `p` at **zero cost**, provided the two
atoms the ripple gives up on are closed through their handles explicitly
(`s10/repairD.py`):

```
a29539 = 12846437·(x_14853 − x_1308) − x_29967   ->  close via x_29967
a40826 (big check, contains x_29967 − x_11360·x_30163)  ->  close via x_30163
```

After that repair the nonzero-atom set is **exactly the same seven** for
`k = 1, 2, 7, −3228258` — nothing else breaks. Since `gcd(p, 7376877) = 1`,
`D mod 7376877` therefore takes every residue, and constraint (1) collapses to

> **(1′)  `A1 + 7376877·A7 ≡ D₀ (mod p)`**

This matters methodologically: earlier sessions' chase reported 29539/40826 as
hard breaks. They are not — the ripple's repair rule (canonical output variable,
exact division only) was simply too weak to find the handle. **Any "the chain is
forced" conclusion reached with that ripple should be re-checked.**

---

## 4. Exact optimum for this defect placement: 5 of 12

With constraints (1′) and (2), enumerate all 2¹² subsets `S` and test integer
solvability exactly — integer kernel of `M_S` by column HNF, then a 2-row integer
linear system for the two congruences (`s10/lattice3.py`):

```
size 12..6 : 0 integer-solvable
size 5     : 300+ integer-solvable
```

> **The optimum is |S| = 5, i.e. 7 failing, i.e. score 39,026 — proved, not searched.**

This confirms Session 9's claim by a cleaner route. Note the earlier *one*-congruence
model (`s10/lattice.py`, treating `x_28730` as free) does admit a 6-subset; the
resulting assignment realises the target atom vector **exactly** but scores 39,011
because `x_28730` breaks 7930/41512. That is a good check that the model is tight.

### The accounting rule

The kernel of `M_S` has dimension `7 − rank(M_S)`, and satisfying `c` independent
mod-p congruences needs `c` free parameters. Hence:

| binding congruences | max |S| | score |
|---|---|---|
| 2 | 5 | **39,026** (current) |
| 1 | 6 | 39,027 |
| 0 | 12 (`A = 0`) | **39,033 — full solve** |

> **There is no partial credit between 39,027 and a complete solution.** One
> congruence is worth exactly one equation; killing both solves the instance.

---

## 5. The MUX branch, measured in equation space

On `x_4287 = 1` (`x_21279 = 1`, `x_7075 = 0`) the two congruences become
satisfiable: `x_2099 = x_37158 + x_9118` and `x_19964 = x_20492 + x_8731` with
both `x_9118`, `x_8731` free. Directly constructed (`s10/muxzero.py`):

> **all seven residual atoms are zero simultaneously** — `a22229 … a35762 = 0`.

The price: 8 collateral atoms (7930, 19088, 22233, 22235, 29539, 37887, 40826,
41512) → **44 failing, score 38,989**.

The MUX's own load pins have handles, but only p-quantised ones
(`s10/dormant.py`):

```
a3568 : x_31861 = C1 + 13479571·x_27676 ,  x_27676 = x_23333·x_6504 = p·x_6504
a3570 : x_14865 = C2 +           x_7574 ,  x_7574  = x_8143 ·x_26658 = p·x_26658
```

so `x_31861 ≡ C1` and `x_14865 ≡ C2 (mod p)` remain pinned. The branch therefore
carries **4 mod-p conditions** (two load obligations `x_4306 ≡ 0`, `x_27177 ≡ 0`
plus the two it just solved) against **2 free residues** — the same deficit of 2,
relocated. This is Session 9 §12.4's conclusion, reproduced here from the
equation-space side and with the load pins' handles accounted for.

---

## 6. Where the remaining freedom is not

- `x_31861`, `x_14865`: move only in multiples of p ⇒ pinned mod p.
- `x_28730`: pinned by `x_4432` through atom 7930.
- `x_7068`: free mod p only (multiples of p), which does not touch `D mod p`.
- The other 9 residual knobs: fully free, and already used optimally.

## 7. Standing next experiments

1. **Re-audit every "forced chain" verdict from Sessions 9 and earlier using
   explicit handle repair rather than the canonical ripple** (§3). This already
   overturned 29539/40826; `s9/chase.py`, `s9/solve_branch.py` and the §14.4
   forced-chain argument all rest on the weaker rule.
2. Atom **7930** (`9367949·(x_24548 − x_25442) − x_7927`) is the single thing
   making congruence (2) binding. Enumerate everything that moves `x_24548`
   and `x_25442` and check whether they can move *together*; if `x_28730`
   becomes free, the score is 39,027 immediately.
3. The two residues `D₀ mod p` and `K2 mod p` remain the whole problem. They are
   `61705020361863629770768910187978745858728889529652486596432934143473517757811`
   and `33310166114805471624282140578459083391052142224394967852279417483154815501175`.

---

## 8. Global handle census — the trapdoor's construction, verified exhaustively

For each atom define its **absorbable lattice**: the subgroup of ℤ generated by
`∂(atom)/∂x_u` over every free input `x_u` whose *only* atom is that one — i.e. the
amounts by which the atom can be shifted at zero cost to everything else.
Computed over all 42,267 atoms (`s10/handles.py`):

```
free inputs occurring in exactly one atom : 1249
  p-quantised, granularity exactly p      : 1240
  dormant handle (rigid)                  :    9
  granularity 1 (unquantised)             :    0
  any other granularity                   :    0
```

> **Every single solo handle in the entire instance has granularity exactly `p`.**
> Not one is unquantised, and none has an intermediate granularity.

That is the trapdoor stated globally and exhaustively, rather than inferred from failed
searches: solo-handle moves can shift any atom only by multiples of `p`, so **every residue
mod p in the system is invariant under them**. In particular neither `D mod p` nor
`K2 mod p` can move, which is precisely why the two congruences of §2 are rigid.
None of the 12 residual equations contains an atom with a free (ℤ) handle.

Session 9 audited ~6 links by hand and reached the same conclusion for those; this is the
same statement for all 1,249.

## 9. Beam search under the *stronger* repair rule

§3 showed the canonical ripple is too weak, so the natural worry is that earlier
"forced chain" verdicts were artefacts. Re-run with the strong rule — effective-linear
solve over **every** variable of a broken atom (`s10/tools.solve_lin`), beam width 200,
depth 10, incremental atom re-evaluation, and the seed variable forbidden as a repair
choice (`s10/beam7930.py`, `s10/beamD.py`).

**The success criterion must be stated carefully.** "Collateral closed" is not enough:
`lib.ripple` recomputes any *gate output* it can, so it will happily restore the seed
variable and report a clean close that changed nothing. The real test is

> collateral empty **and** the target residue (`D mod p` or `K2 mod p`) actually moved.

Results:

| seed | collateral closes? | residue moved? |
|---|---|---|
| `x_28730 += k·p` (k = 1, 2, −5) | **yes**, depth 2 via `x_7927` then `x_11052` | no — `K2` unchanged, `a22230` moves by exactly `k·p` |
| `x_28730 += 1, −1, 2, 7376877, p+1, 12846437` | no — ladder 11625 → 11624 → 11621 → 30238 → 24948 → … | — |
| `x_7068 += ±1, +12846437`; `x_7068 +1` with `x_14853 +1`; `x_10878 += 1` | no — ladder 27314 → 29539/40826 → 18686/39719 → 19482 → 19480 → … | — |
| `x_2099 += 1`, `x_37158 += 1`, `x_22542 += 1` | "yes" | **no — artefact**: the ripple recomputed `x_2099` from its definer 29090 and put it back. `D mod p` identical, score identical to base (39,026) |

So the Session-9 conclusion **survives the stronger repair rule**: every move that closes
is `≡ 0 (mod p)` and leaves both residues fixed. What §3 corrects is not the verdict but
the *evidence* — the p-multiple moves genuinely do close and the old ripple said they did
not, so any argument that relied on those breaks being real needs redoing.

## 10. Summary of the session

* Deliverable re-verified and unchanged: **39,026 / 39,033**.
* The residual is now in exact closed form: 7 atoms, 12 equations, rank 7, and an
  achievable set described by exactly two mod-p congruences.
* **39,026 is proved optimal for this defect placement** by exhaustive exact integer
  enumeration — not by search.
* The accounting rule says one congruence is worth exactly one equation and killing both
  is a complete solve: **39,026 → 39,027 → 39,033**, nothing in between.
* The global handle census shows every one of the 1,249 solo handles is exactly
  p-quantised, which is the structural reason both congruences are rigid.

---

# Part II — the global attack (same session)

Deliverable still **39,026**. This part abandons local repair entirely and attacks
the instance globally. It produces the first structural crack in the design.

## 11. The forward-eval frame — the honest global statement

Take the 39,026 witness's FREE INPUT values and forward-evaluate every gate in
topological order (`s10/forward.py`). Result:

> **6 nonzero atoms, all CHECKS, zero broken gates** — 37 failing equations (38,996).

```
a7930   9367949*(x_24548 - x_25442) - x_7927     ->  x_24548 == x_25442 (mod p)
a29539  12846437*(x_14853 - x_1308) - x_29967    ->  x_14853 == x_1308  (mod p)
a35759  5113045*x_7075*x_9118 - x_29854          ->  x_9118  == 0       (mod p)
a35760  x_31864 - x_28961*x_10903                ->  x_8731  == 0       (mod p)
a40826, a41512                                    big checks, 1 equation each
```

Every one contains a **free input** (`x_24548`, `x_14853`, `x_9118`, `x_8731`).
This also explains the 39,026 witness: it deliberately violates five GATE atoms
(22229, 22230, 35758, 35761, 35762) precisely so that `x_1308` and `x_25442` land
on `x_14853` and `x_24548`. The seven-atom picture of Part I is that trade.

## 12. Exact machinery: reverse-mode AD mod p

`s10/ad.py` propagates adjoints backwards through the gate DAG
(`dx_t = -sum_w (da/dx_w)/(da/dx_t) dx_w`), giving `d(check)/d(free input)` mod p
in one pass per check. **Validated against exact finite differences** — every
non-boolean input matches; only the 0/1 controls mismatch, as they must.

Gradient supports are tiny: 2, 5, 9, 80, 132 free inputs.

## 13. The point is RIGID — measured, not inferred

A step must zero the failing checks *and* keep the ~10,786 satisfied ones. Only
checks reachable from the step's support can move, so the correct system closes
(`s10/closure.py`):

```
round 0: cols=12  rows=50   round 3: cols=75 rows=187
round 4: cols=79  rows=193  -> FIXED POINT
```

`s10/rankdef.py`: **rank(A) = 79 of 79 columns — full column rank, zero null
space — with 6 independent inconsistencies.** Degenerate (square-check) rows: 0,
so this is not Session 9's Jacobian artefact.

Relaxing all 256 message bits to arbitrary GF(p) values closes at
**2,352 rows x 710 columns (exactly 256 boolean) and is still inconsistent**
(`s10/closure_bits.py`). A full single-flip scan of all 1,156 boolean free inputs
with genuine forward evaluation (`s10/bitscan.py`) finds nothing better than 37.

> Local and first-order methods are definitively dead. Any solution requires a
> structural change.

## 14. THE CRACK: the p-wire is not rigid

Every handle in the instance enters as `wire * handle`, where `wire` is one of
**220 variables all equal to p**. That single fact is why the Part I census found
all 1,249 handles p-quantised. The wire is 219 copies of one root, held by a
single **bare pin** `a37694 = x_26064 - p`, which appears in only **12 equations**
(`s10/wire.py`).

Setting the whole wire to 1 (`s10/wire1.py`) flips the census completely:

| | granularity p | granularity 1 |
|---|---|---|
| p-wire | **1240** | 0 |
| wire = 1 | 0 | **1240** |

On that branch the congruences dissolve — `a7930` and `a29539` close exactly
through their handles `x_11052`, `x_30163`, and the two big checks `a40826`,
`a41512` come along for free (`s10/trade2.py`):

> **wire = 1 reaches 39,020 with only TWO nonzero atoms** — the wire pin `a37694`
> (12 equations) and `a39417` (1 equation). Its 13 failing equations contain
> **only wire-copy atoms and boolean pins — no free inputs at all.**

## 15. The wire deformation kernel — dimension 3, not 0

Write `w_u = p + d_u`. Every wire-identity atom is then linear and homogeneous in
`d` (copy atom `x_i - x_j` -> `d_i - d_j`; root pin -> `d_root`), so all 219
equations containing them become a homogeneous system `M d = 0` in `Z^220`
(`s10/wirekernel.py`):

```
rank(M) = 217 of 220     ->     KERNEL DIMENSION 3
```

> **The wire is NOT rigid. There are three directions in which it can deform
> without breaking a single wire-identity equation.**

Per-coordinate reachability (`gcd` of the basis at each coordinate):

```
members the kernel can move        : 217 of 220
members with gcd 1 (ANY value)     : 161
  including handle multipliers x_11360, x_28599, x_17499, x_22665, x_28961
x_15616 : gcd 29        x_26064 (the root) : gcd 0  -- FIXED
```

So five of the six handle multipliers can be set to **any value, including 1, for
free** — which would unquantise their handles outright. The root pin is the one
thing the kernel cannot touch.

## 16. Why the deformation does not (yet) pay

Applying a kernel vector and re-solving every handle to restore the original gate
outputs (`s10/deform2.py`) succeeds for **3,346 of 3,349** product gates. The
damage classifies as (`s10` diagnostic):

```
235 broken atoms
  215  wire copy atoms  <- EXPECTED; their equations still cancel by construction
   13  multi-wire monomials (w_i * w_j)
    3  downstream, 3 single-wire gates, 1 check
```

so the genuine cost is **20 atoms**, matching Session 8's "~13 active unpackings".
Net measured score on that branch: 38,981, and about 39,018 once the six checks
are then closed — still short of 39,026. The obstruction has moved from the
handles to the **multi-wire monomials `w_i·w_j`**, whose invariance is *quadratic*
in `d` (`p(d_i + d_j) + d_i d_j = 0`) and therefore not captured by the linear
kernel.

## 17. Where this leaves the attack

The design's p-quantisation — the thing Part I proved makes both congruences
rigid — **is breakable**: the wire has a 3-dimensional free deformation space and
161 of its members can take any value at zero cost to the wire-identity
equations. What now stands in the way is a much smaller and sharper object:

> **13 multi-wire monomials `w_i·w_j`, and the single bare root pin `a37694`
> (12 equations) that the kernel cannot move.**

Highest-EV next experiments, in order:
1. Solve the deformation with the multi-wire monomials imposed **exactly**
   (`p(d_i + d_j) + d_i d_j = 0` — a quadratic system in 3 unknowns after
   restricting to the kernel). If it has a nonzero solution, the handles
   unquantise at zero cost and both congruences fall.
2. Equation-space compensation for the 13 equations of `a37694` using deformed
   copy atoms — they are all wire-identity atoms, so this is the same linear
   algebra with the pin's row moved to the right-hand side.
3. LLL-reduce the 3-dimensional kernel lattice: the current basis has ~325-digit
   entries; a short vector would make the whole branch numerically tractable.

## 18. The deformation route, closed — and exactly why

`s10/multiwire.py` resolves §16's remaining obstruction. The multi-wire monomials
are almost all **self-pairs** `w_i·w_i`: they come from the degree-4 **square**
check atoms, where a wire variable appears squared. Their invariance condition
`w_i·w_j = p²` therefore degenerates to

> **`w_i² = p²`  ⇒  `w_i = ±p`.**

The constraint graph over these pairs is a **single component spanning all 220 wire
members, non-bipartite, with self-loops**. So every wire member is forced back to
`±p`, and `−p` leaves the handle granularity at `p` — no gain.

> **The wire-deformation crack is real but bounded: the linear wire-identity
> equations permit a 3-dimensional deformation, and the ~20 degree-4 square checks
> then pin every member back to ±p.** Each of those checks lives in exactly ONE
> equation, which is why the branch measures 38,981 (≈39,018 after closing the six
> checks) rather than collapsing entirely — the price is ~20 single-equation checks,
> and that price exceeds the 7 equations the p-wire branch already pays.

### Revised standing recommendation
The obstruction is now located to the sharpest object found in any session:

* the **linear** wire-identity system has a 3-dimensional kernel (the wire is not rigid);
* the **quadratic** square checks `w_i² = p²` are what actually pin it;
* the root pin `a37694` (12 equations) is untouchable by the kernel.

Any further attack should aim at the ~20 degree-4 square checks that carry a squared
wire variable — they are individually cheap (1 equation each) and are the only thing
standing between the 3-dimensional deformation and an unquantised handle set. A
deformation that keeps `|w_i| = p` for exactly the members appearing squared, while
moving the members that only ever appear linearly, would cost nothing at all; whether
the handle multipliers `x_11360, x_28599, x_17499, x_22665, x_28961` are among the
squared ones is the single question to settle first.

### 18.1 The question settled — the route is comprehensively closed

Counting, over all 42,267 atoms, which wire members occur with multiplicity ≥ 2 in
some monomial:

```
wire members appearing SQUARED : 220 of 220
wire members NEVER squared     : 0
  x_11360  squared in 3 monomials      x_22665  squared in 8
  x_28599  squared in 9                x_28961  squared in 9
  x_17499  squared in 9                x_15616  squared in 3
```

> **Every one of the 220 wire members is squared somewhere, so every one is pinned to
> `±p`. There is no subset that can move for free.** The setter covered the whole wire
> with degree-4 square checks; that — not the linear identity chain — is what makes the
> p-quantisation rigid.

This is the correct final statement of the design, and it supersedes §17's optimism:
the 3-dimensional linear kernel is real, but it is entirely absorbed by the quadratic
square-check constraints. The p-quantisation of all 1,249 handles therefore stands,
and with it the two mod-p congruences of Part I and the optimality of **39,026**.

---

# Part III — the budget attack, and the trapdoor priced exactly

Deliverable still **39,026**. This part reframes the instance adversarially and
produces the sharpest characterisation of the design in any session.

## 19. The reframe: a 6-equation budget

The current branch pays **7** failing equations, so *any* structural violation
costing ≤ 6 beats it. That turns every guard into a price tag
(`s10/pricelist.py`, price = how many equations a check lives in):

| guard | price |
|---|---|
| a degree-4 square check `E²` | **1 equation** |
| `a40826`, `a41512` | **1 each** |
| a wire copy atom | 12–14 |
| the wire root pin `a37694` | 12 |
| a boolean check `x²−x` | 13–15 |

Two mistakes in my earlier greedy searches, now fixed: they scored by *number of
nonzero atoms* rather than failing equations, and refused any move that raised the
atom count even when the new atoms were 1-equation checks and the closed one cost
15. Re-run with the correct objective (`s10/pricelist.py`, beam over equations):
best 39,005 — still short, but the objective is now right.

## 20. The six inconsistency certificates

Augmenting the closed system as `[A | b | I]` and eliminating on `A` extracts
explicit left-null vectors `y` with `y·A = 0`, `y·b ≠ 0`, and names the checks
that combine to produce each (`s10/certs.py`):

```
rank(A) = 79 of 79        INCONSISTENCY CERTIFICATES: 6

cert 0: 20 checks; cheapest members (1, 36602) (1, 37887) (7, 35759) (8, 35760)
cert 1: 12 checks; cheapest members (10, 2423) (10, 21617) (10, 31670) (11, 19297)
cert 2: 13 checks; cheapest (1, 41400) (10, 2423) (10, 31670) (10, 34397)
cert 3: 11 checks; cheapest (1, 41507) (10, 2423) (10, 31670) (10, 40065)
cert 4: 13 checks; cheapest (1, 41827) (10, 2423) (10, 31670) (10, 34397)
cert 5: 19 checks; cheapest (1, 11007) (1, 25676) (1, 39800) (1, 42245)
```

> **Five of the six certificates can be hit for 1 equation. Certificate 1 cannot —
> its cheapest member costs 10, and its members `2423, 21617, 31670, 19297` are
> exactly the §15.2 trapdoor chain.**

Minimum-cost hitting set (greedy + swap): `{21617, 36602, 41400, 41507, 41827,
42245}` at **15 equations**. Note `2423`, `31670` and `19297` each hit *four* of
the six certificates, so the real shape of the optimum is `cost(hub) + 1 + 1`.

> **The design carries a margin of exactly 8: the cheapest way through costs 15,
> the give-up option costs 7.**

## 21. Hub compensation does not pay — and why (a corrected claim)

`s10/truecost.py` first suggested hub `a31670` had true cost **1** (nine of its ten
equations appeared to have a compensating atom). **That heuristic was wrong**: it
counted the *existence* of an adjustable helper per equation, not the fact that a
helper's value is a single number shared across all of them. The exact computation
(`s10/hub31670.py`) settles it:

```
a31670 = (x_22152 - HUGE) - 7550763*x_29309       a31669 = x_29309 - p*x_105
=> (a31670, a31669) = (D - 7550763*s, s - p*h)  -- a31669 free, a31670 fixed mod 7550763
region matrix over its 10 equations + helpers: 41 equations x 16 atoms, RANK 16 (full)
```

With `a31670 ≠ 0` forced there are ~2 free parameters against 10 equations, so at
most **one** can be saved; and recruiting the other 15 adjustable atoms drags **31
additional equations** into the region, costing far more than it saves. Hub cost is
therefore ≈ 9, not 1, and the 15-equation hitting set stands.

## 22. Other lines priced and closed this round

* **Wire deformation with square-check repair** (`s10/deform3.py`): my §18 closure
  used a *sufficient* condition (monomial invariance `w_i·w_j = p²`), not a
  necessary one, so I re-tested whether the free variables inside each square check
  `E` can absorb a wire change. Across all six kernel directions the repair found
  **zero** admissible moves — after deformation the handles carry ~325-digit values
  and the exact-division condition fails. §18's conclusion survives the stronger test.
* **The forced OR gate is a non-lever.** Every load pin `bit·(x_B − HUGE)` is free
  to satisfy when its bit is 0, so the HUGE constants — the only source of values
  that are not multiples of p — enter *only* through set bits. But `x_9274 = 1`
  holds automatically at all-bits-zero, and that branch measures **38,871**
  (31 nonzero checks). Only two boolean free inputs are set at the deliverable
  (`x_2081`, `x_24601`), so the instance is already almost all-zero on the bits.

## 23. Standing assessment

The instance is now priced rather than merely characterised:

* the obstruction is exactly **6 independent certificates**;
* **5 of them are cheap** (1 equation each — the single-equation square checks);
* **certificate 1 is the trapdoor**, minimum price 10, and it is the same chain
  §15.2 audited link by link;
* the total is **15 against a budget of 7** — an 8-equation margin.

The productive question is no longer "can the system be solved" but **"can
certificate 1 be hit for under 9 equations?"** Everything else is already cheap.
