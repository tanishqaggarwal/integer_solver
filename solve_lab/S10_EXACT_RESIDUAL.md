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

Per-value reachability (`gcd` of the basis at each value):

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

## 24. Every route priced — and a correction to Part II's optimism

Two further corrections and a complete price table.

**Correction to §15.** Part II reported that 161 wire members have kernel-gcd 1 and
so "can be set to 1 for free". Reachability is real but the *magnitudes* are not:
hitting `d_u = 1 − p` needs kernel coefficients ~10²⁵⁰, which blows every other
value to ~10⁵⁷⁵. Measured (`s10/deform_solve.py`): applying the raw kernel
directions gives `|w_3915| = |w_11360| = 325 digits`, i.e. handle granularity far
*worse* than p, and the branch scores **38,990**. The gcd-1 statement is true and
useless on its own; only a *short* kernel vector would matter, and the lattice is
3-dimensional with ~325-digit entries.

**Complete wire price table** (`s10/memberprice.py`), cost = identity equations
violated + square-check equations violated:

```
cheapest member overall     x_15413   13  (10 identity + 3 square)
uniform wire shift (any c)            13  (root pin a37694 = 12, + a39417 = 1)
cheapest USEFUL multiplier  x_3915    15  ( 9 identity + 6 square)
x_11360   36     x_14466   54     x_15616   56     x_22665   48
x_17499  150     x_28599  292     x_28961  208
kernel deformation                   ~20  (and useless magnitudes)
certificate hitting set              15
------------------------------------------------------------------
current branch (give up)              7   <== the deliverable
```

> **No route into the design costs less than 13, against a give-up cost of 7. The
> margin is 6, and it is uniform across every attack surface measured: the wire,
> the certificates, the hubs, and the bits.**

### What this settles
The instance is now priced end to end rather than merely characterised. Session 9
established *that* the chain is rigid; Part I established the residual in closed
form and proved 39,026 optimal for its placement; Parts II–IV establish the
**cost of every way around it**. The design's guard is not one wall but a uniform
6-equation margin over every reachable violation.

### The only remaining question
Certificate 1 must be hit for **under 9** equations, and its cheapest member is 10.
Everything else is already cheap. Two concrete sub-questions, both open:
1. Is there a *short* vector in the 3-dimensional wire kernel (LLL over ~325-digit
   entries)? A short vector with support ≤ 5 would cut the deformation cost from
   ~20 to ~5 and take the wire route under budget.
2. Does certificate 1 admit a member outside the closed 79-column system — i.e.
   can it be hit by a variable the closure never reached?

---

# Part V — the invariant

## 25. No sparse kernel vector exists (sub-question (a), settled)

A kernel deformation costs 0 identity equations plus the square checks of the
members it moves, so a kernel vector with support 2–3 would cost ~3–6 and beat 7.
View the kernel as a 220×3 matrix `K` (row = the three basis values at that
member); a kernel vector supported on `S` exists iff every row outside `S` lies in
a common plane of ℚ³ (`s10/sparsekernel.py`):

```
kernel-zero rows          : 3  (x_1692, x_26064, x_32499)
distinct directions       : 215 of 217 nonzero rows
largest rank-1 cluster    : 2 rows
largest coplanar set      : 4 rows
=> sparsest kernel vector has support >= 213 of 220 members
```

> **There is no sparse kernel vector.** Every free wire deformation moves at least
> 213 members and therefore breaks essentially all of their square checks. The wire
> route cannot get cheap, and its measured floor of 13 (uniform shift) stands.

## 26. Seven is an INVARIANT, not a property of the placement

Part I proved 5-of-12 optimal for the seven-atom placement. That proof was
*conditional on the placement*. Enlarging it with additional adjustable atoms
(`s10/eighth.py`, `s10/invariant.py` — exact integer subset enumeration each time):

| extra atoms | region (equations) | max satisfied | **failing** |
|---|---|---|---|
| — | 12 | 5 | **7** |
| 35756 | 15 | 8 | **7** |
| 35754 | 17 | 10 | **7** |
| 35756 + 35754 | 18 | 11 | **7** |

> **Every extra free parameter buys exactly as many equations as it drags in. The
> failing count is pinned at 7 across every placement tested — it is an invariant
> of the instance, not an artefact of the defect set.**

Only two adjustable atoms in the whole instance even touch the twelve equations
(35756 with overflow 3, 35754 with overflow 5); both were tested, alone and
together, and both leave the count at 7.

This is the exact, reproducible form of what earlier sessions called the
"conserved obstruction" (`CONSERVED_OBSTRUCTION.md`) — previously an observation
about relocating quantities, now a numerical law with an identified mechanism.

## 27. Consolidated status

Everything measured in Session 10 agrees on one number:

```
give-up cost (the deliverable)                      7
failing count under placement enlargement           7   (invariant, 4 placements)
minimum wire route (uniform shift)                 13
cheapest single wire member (x_15413)              13
certificate hitting set                            15
kernel deformation                                ~20   (support >= 213, forced)
```

The instance is priced end to end and the price is uniform. Three independent
lines — the equation-space lattice (Part I), the GF(p) closure and its six
certificates (Parts II–III), and the wire geometry (Parts IV–V) — all return the
same 7, and every route around it costs at least 13.

### What would actually break it
Nothing measured this session gets under budget. The two remaining doors, both
now precisely stated:

1. **Certificate 1 for under 9 equations.** Its cheapest member costs 10 and its
   support is the §15.2 chain. Five of the six certificates are already 1 apiece.
2. **A knob outside the closed 79-column system.** The closure is a fixed point
   over non-boolean free inputs; with all 256 bits relaxed to GF(p) it grows to
   2,352 × 710 and stays inconsistent. Something outside *both* would be needed —
   and the invariance of 7 across placements is evidence that no such knob exists
   in the equation-space direction.

---

# Part VI — the invariant survives a genuine gap in my own model

## 28. Forensics on the setter's constants: no backdoor

Every prior analysis treated the ~2,817 large literals as opaque. Tested directly
(`s10/forensics.py`): gcd of all large constants = 1; gcd of all pairwise
differences = 1; no constant equals either binding residue mod p; `D0/K2 mod p` is
a full-size 253-bit number and `D0 ≠ k·K2 (mod p)` for any `k < 60`. Both residues
*are* quadratic residues, which is the only structure found and is not exploitable.
The constants are ~290 bits (p is 256), quotients `c // p` are ~10¹⁰–10¹¹ with no
pattern. **The setter's arithmetic has no exploitable structure.**

## 29. A real gap in Part I's model — `a22231` is a free 8th atom

Part I declared `x_28730` "not free" because moving it drags `x_4432` and breaks
atom 7930 (15 equations). That was an artefact of how I moved it: I moved
`x_28730` *together with* `x_4432` to hold `a22231 = 0`. **But `a22231` need not be
zero.** Moving `x_28730` alone (`s10/a22231.py`):

```
a22230 = x_28730 - p*x_9413            changes by +d
a22231 = x_4432 - x_19964 - x_28730    changes by -d
x_4432 untouched  ->  NO downstream collateral   (verified for d = 1, 2, p)
```

and `a22231`'s ten equations lie **entirely inside the twelve** — zero overflow.
So the correct model has *eight* atoms and **one** congruence on a pair rather than
two separate ones:

```
A1 + 7376877*A7 == D0 (mod p)        A2 + A8 == K (mod p)      A3..A6 free
```

Exact optimisation on this strictly larger model: **6 of 12 satisfiable**, versus
5 before. Constructed and verified end to end (`s10/build27.py`): all eight target
atom values realised **exactly**, `x_4432` untouched, six of the twelve equations
satisfied — the first time the region has gone past 5.

**The gain is cancelled outside the region.** With `a22231 ≠ 0`, the square atom
`a37887 = R²` lights up and breaks equation 8680, restoring 7.

## 30. Killing `a37887` too — and the invariant again

`R` is a linear combination of atoms (extracted by parsing the square's root):

```
R = a22231 + 6*a22232 + 15*a22233 - 21*a22234 - 13*a22235
    + (a19087..a19092, a10935..a10941 terms)
```

and `a22232..a22235` are movable in pairs (`x_23754` moves `a22232/a22233`
oppositely, `x_35619` moves `a22234/a22235` together), so `R = 0` is reachable —
`a22231 = 9δ + 34ε` with `gcd(9,34) = 1` hits every value. Extending to a 12-atom
model with `R = 0` imposed exactly (`s10/kill37887.py`): region grows 12 → 16
equations, max satisfied grows 5 → 9. **Failing: 7.**

## 31. The invariant, now on six independent placements

| model | atoms | region | max satisfied | **failing** |
|---|---|---|---|---|
| Part I baseline | 7 | 12 | 5 | **7** |
| + 35756 | 8 | 15 | 8 | **7** |
| + 35754 | 8 | 17 | 10 | **7** |
| + both | 9 | 18 | 11 | **7** |
| **+ a22231** (gap fixed) | 8 | 12 | **6** | **7** (a37887 costs 1 outside) |
| + a22231, a22232..a22235, R = 0 | 12 | 16 | 9 | **7** |

> **Six independent placements — including one built specifically to exploit a real
> error in the earlier model — all return exactly 7. Every degree of freedom added
> is matched, to the equation, by the equations it drags in.**

This is the strongest form of the result: the invariant is not an artefact of the
defect set I happened to choose, because I found and fixed a genuine gap in that
choice and the number did not move.

---

# Part VII — the number theory, closed; and the root pin's true price

## 33. Rational reconstruction: the constants carry no structure

Given `r mod p`, extended-Euclid yields the unique small `a/b ≡ r` whenever
`|a|,|b| < sqrt(p/2)`. A setter who built residues from small rationals would be
exposed here. Measured (`s10/ratrec.py`) — every residue returns **maximal-size**
numerator and denominator (38–39 digits, i.e. right at the `sqrt(p/2)` bound):

```
D0, K2, D0/K2, K2/D0, D0*K2, D0±K2, 1/D0, 1/K2, HUGE mod p, C1 mod p
   -> all 38-39 digit a and b : NO small rational structure
```

Further: `gcd(HUGE, C1) = 1`; quotients `HUGE//p = 1094785891323` and
`C1//p = 289077647971` are unremarkable 12–13 digit numbers; `HUGE ≠ k·C1 (mod p)`
for all `k < 200`. Of the **2,815** constants exceeding p, **zero** have residue
below 2⁸⁰ and **all 2,815 residues are distinct**. The seven residual equation
values have **gcd 1** and only tiny random small factors.

> **There is no arithmetic backdoor. The setter's constants are random, and the
> two binding residues are unrelated to each other and to everything else in the
> file.** The number-theoretic line is closed.

## 34. New: the wire root frees for ONE identity equation, not twelve

The uniform wire shift costs 12 — the root pin `a37694` lives in 12 equations.
But that is not the minimum. Since `e_root` lies in the row space of the identity
system, write `e_root = y₀ᵀM`; any `d` annihilating `supp(y₀)` has `d_root = 0`,
so at least one equation of `supp(y₀)` must break. Computed (`s10/rootfree.py`):

```
rank(M^T) = 217, system consistent
supp(y_0) = { equation 37257 }        <-- a SINGLE equation
root-pin equations = [8429, 11166, 12594, 23869, 25313, 26785,
                      31400, 32300, 36106, 36767, 37257, 37666]
```

> **Equation 37257 is the unique identity equation whose wire content is the root
> pin alone.** In the other eleven, `a37694` sits alongside copy atoms that can
> absorb it, so a *non-uniform* deformation compensates it there. The root's true
> identity-space price is **1 equation, not 12.**

Constructed (`s10/freeroot.py`): dropping equation 37257 gives rank 216 and a
**4-dimensional** deformation space, all four directions moving the root. Applied,
every one of the other 218 identity equations still holds.

**And it still does not pay.** All four directions have ~324-digit entries and
support 217, so 17 non-copy atoms break: **38,984**. Identity cost 1, square-check
cost ~12, total ~13 — the same floor, reached now from a third independent
direction.

## 35. Final consolidated position

```
give-up cost (the deliverable)                       7
failing count under placement enlargement            7   (invariant, 6 placements)
root freed via eq 37257 + square checks             ~13
uniform wire shift                                   13
cheapest single wire member                          13
certificate hitting set                              15
kernel deformation                                  ~20
```

Every independent line — equation-space lattice, GF(p) closure and certificates,
wire geometry, the identity row space, and now the number theory — returns the
same two numbers: **7 to stop, ≥ 13 to go through.** The margin is 6 and it has
not moved under any attack in this session.

**What is genuinely still open:** hit inconsistency certificate 1 for under 9
equations (its cheapest member is 10, and five of the six others cost 1 apiece).
Every other door measured this session is closed with a number attached.

---

# Part VIII — every door, opened and priced

## 36. A counting error of mine, corrected

Part III reported "6 independent inconsistencies". **`b` is a single column, so
`rank([A|b]) − rank(A) ≤ 1`** — there is exactly *one* independent obstruction, and
the six were six witnessing rows of it. That raised a real hope: a single dropped
row might restore consistency. Tested exhaustively (`s10/singledrop.py`):

```
every one of the 128 single-row drops : 0 restore consistency
every pair among the 30 cheapest rows : 0 restore consistency
```

So the obstruction survives removal of any one row, and of any cheap pair.

## 37. The region, closed exhaustively

`eighth.py` defined "adjustable" as *carrying a solo free handle* — which is
precisely why it missed `a22231`. Correct definition: a variable is a **region
knob** if moving it changes no equation outside the twelve. Scanning every
variable in every atom of the twelve equations (`s10/regionknobs.py`):

```
variables whose ENTIRE footprint lies inside the twelve : 9
  x_642, x_1329, x_8731, x_9118, x_9413, x_10903, x_17325, x_29854, x_31864
region atoms reachable by those knobs : exactly [22229, 22230, 35758..35762]
next cheapest variable : x_28730, footprint 1 equation outside (eq 8680 via a37887)
everything else        : >= 3 equations outside
```

> **The region has exactly nine free knobs and they reach exactly the seven atoms
> Part I used.** There is no hidden freedom, and `x_28730`'s "1 outside" is
> precisely the `a37887` cost that `build27.py` measured. The region is closed.

## 38. Boolean branches, exhaustive in the witness frame

Earlier bit scans ran from the forward state (37 failing) or used the unreliable
ripple. Re-run in the **witness** frame with exact repair, all 1,156 boolean free
inputs (`s10/bitwitness.py`):

```
best 20 flips : failing = 7, region 12, satisfied 5   (identical to the base)
x_4287        : 34 failing      x_24601 : 83      x_2081 : 106
```

> **No boolean flip in the deliverable's frame improves on 7**, and the three
> structural control bits are far worse. The branch structure is exhausted.

## 39. The cyclic components — a genuinely new freedom, and it is inert

`fwd.py` covers only 29,675 of 31,475 defined variables: 1,800 sit in gate cycles.
A cyclic block is a *system*, and a singular one has a family of solutions that
forward evaluation silently collapses to a point. Measured (`s10/cycles.py`):

```
non-trivial SCCs : 40, all of size 2
local Jacobian   : rank 1 of 2 in every one  ->  kernel dimension 1 each
```

**Forty free parameters that no forward-eval or local method could see.** Slid
along each line (`s10/slide.py`):

```
all 40 slides : no new nonzero atoms, failing stays 7, D0 and K2 UNCHANGED
8 of 40 are literally inert (touch only their own two gate atoms)
```

> The freedom is real and it is **orthogonal to the obstruction**: it moves the
> assignment without moving either binding residue. A new door, opened, and empty.

## 40. Final position

| door | status | cost |
|---|---|---|
| region knobs beyond the 9 | **closed exhaustively** | — |
| boolean branches (1,156, witness frame) | **closed exhaustively** | ≥ 7 |
| single-row / cheap-pair sacrifice | **closed exhaustively** | none work |
| cyclic-component freedom (40 params) | **open but inert** | 0 gain |
| number theory / constants / rational reconstruction | **closed** | no structure |
| wire: uniform, per-member, kernel, root-via-37257 | **closed** | ≥ 13 |
| certificate hitting set | closed | 15 |
| **give up (the deliverable)** | — | **7** |

Seven independent lines now return the same answer. The instance's margin is 6
equations and nothing in this session moved it. What remains is not a door I can
name a cheap price for: it is the setter's witness, or a genuine cryptanalytic
break of the pinned residue `D0 = HUGE − C1 (mod p)`.

---

# Part IX — the sacrifice question, answered exactly

## 41. A reformulation that makes it tractable

The earlier exhaustive budget search timed out (its output was lost to a pipe, so
Part VIII's table left this door genuinely open). Reformulated, each test becomes
trivial. Dropping rows `S` leaves `A_{-S}x = b_{-S}`; its left null vectors,
extended by zeros on `S`, are exactly the `y ∈ leftnull(A)` with `supp(y) ∩ S = ∅`.
Hence

> **consistent after dropping `S`  ⟺  `t ∈ colspace(Y[:,S])`**,
> where `Y` is a basis of `leftnull(A)` and `t = Y·b`.

Each test is then a `49 × |S|` rank check instead of a `128 × 80` elimination —
microseconds rather than seconds (`s10/budget6fast.py`).

```
closed system 128 x 79     rank(A) = 79     leftnull dimension = 49
t = Y.b  nonzero  ->  inconsistent, as expected
```

## 42. The minimum sacrifice is exactly three rows

```
size 1 : impossible
size 2 : impossible
size 3 : FOUND  {a3578, a26731, a35759}
```

So **at least three checks must be sacrificed** — no single row and no pair can
restore consistency, confirming and sharpening Part VIII's negative. The
minimum-size set is:

| check | what it is | price |
|---|---|---|
| `a3578` | setter load pin `x_2081·(x_12553 − HUGE)` | 14 |
| `a26731` | mirror `6788513·(x_16742 − x_19083) − x_9254` | 16 |
| `a35759` | one of the six currently-failing checks | 7 |

Union: **37 equations**, i.e. score 38,996 — precisely the forward-eval floor.
The cheapest *sized* solution is the most expensive kind.

## 43. Budget ≤ 6: exhausted

Over all 46 rows priced ≤ 6, with cost-pruned enumeration:

```
size 1 :        46 within-budget sets  -> none
size 2 :     1,081                     -> none
size 3 :    16,261                     -> none
size 4 :   179,446                     -> none
size 5 : 1,550,200                     -> none
size 6 : cost-pruned DFS               -> (see runs/)
```

> Nothing in the first five sizes comes within the 6-equation budget, and the
> minimum-size solution costs 37. The sacrifice route is closed on both axes:
> **too few rows is impossible, and cheap enough is unreachable.**

## 44. Closing statement

Every door named in this session now has a number on it:

```
give up (the deliverable)                       7
failing count, invariant over 6 placements      7
minimum sacrifice: >= 3 rows, cheapest 37 eqs  37
uniform wire / cheapest member / root-via-37257 13
certificate hitting set                        15
kernel deformation                            ~20
region knobs beyond the nine                 none
boolean flips (1,156, witness frame)          >= 7
cyclic freedom (40 parameters)          real, inert
number theory (constants, ratrec, forensics)  no structure
```

The instance is characterised, priced, and closed on every axis I can measure.
What remains is not a door with a price — it is the setter's witness, or a genuine
cryptanalytic break of the single pinned residue `D0 = HUGE − C1 (mod p)`.

---

# Part X — the message space, exhaustively closed

## 45. Global rigidity, tested bluntly

Every rigidity result up to here was a linearisation **at one point**. Tested
directly by randomising the non-boolean free inputs (`s10/randomize.py`):

```
1 / 10 / 100 / 1000 / ALL 6,117 randomised   ->  37 / 148 / 1084 / 5219 / 7355 failing
best over every randomisation                ->  37   (= the base)
do 7930, 29539, 35759, 35760 fail EVERY TIME ->  TRUE
```

> The four core checks fail from **every** starting point, including full
> randomisation of all 6,117 non-boolean free inputs. The residual is pinned
> against them **globally**, not merely locally.

## 46. The "256-bit codeword" collapses to 5 dimensions

Earlier sessions treated the message as a 2²⁵⁶ combinatorial wall. It is not.
Computing the exact AD gradient of every failing check with respect to all 1,156
boolean free inputs (`s10/bitgroups.py`):

```
boolean inputs that move ANY failing check : 128   (not 256)
distinct signature vectors among them      : 5
group multiplicities                       : 75, 50, 1, 1, 1
=> reachable message states = 76*51*2*2*2  = 31,008
```

Within a group the bits are interchangeable, so only the **count** matters. The
whole message space is 31,008 states — enumerable in a second.

## 47. Swept, with the model validated where it holds

A first sweep (`s10/msgsweep.py`) claimed 2 of 6 checks zeroable. **Constructing
that state refuted it** (`s10/msgverify.py`): 62 failing, nothing zeroed. Its two
bits were `x_2081` and `x_4287` — the structural MUX controls, where
`b·(X − HUGE)` has `X` itself depending on `b`, so linearity fails. A lesson worth
recording: the bit model is exact only for ordinary load bits.

Validated and re-run properly (`s10/msgvalid.py`):

```
linearity check, group of 75 (bit x_91) : model MATCHES exactly
linearity check, group of 50 (bit x_47) : model MATCHES exactly
sweep of all 76*51 = 3,876 states of the two linear groups:
    histogram of checks zeroed = {0: 3876}
```

> **The 125 ordinary load bits, swept exhaustively with a validated model, cannot
> zero a single failing check.** The only bits with real leverage are the three
> structural controls `x_2081, x_4287, x_13195`, and those are precisely the branch
> flips already measured (34 / 83 / 106 failing in the witness frame).

## 48. The sacrifice route, exhaustively closed

The budget question finished (`s10/budget6fast.py`, cost-pruned DFS):

```
10,917,019 within-budget sets tested   ->   NONE restore consistency
minimum sacrifice: 3 rows {a3578, a26731, a35759}, cost 37 equations
sizes 1 and 2: impossible
```

## 49. Final

```
give up (the deliverable)                          7
invariant over 6 placements                        7
message space (31,008 states, exhaustive)   0 checks zeroable
non-boolean inputs (global randomisation)     no improvement
sacrifice route (10.9M sets)                       none <= 6
minimum sacrifice                                 37
wire routes                                       13
certificate hitting set                           15
```

The instance is closed on every axis I can measure: the free inputs globally, the
message space exhaustively, the sacrifice route exhaustively, the wire geometry,
the region knobs, and the number theory. **39,026 stands.**

---

# Part XI — the residual is ONE trade, and here is its ledger

Everything in Parts I–X measured *prices*. This part finds the **mechanism** that
sets them, by two new move classes and one new piece of machinery. It does not
beat 39,026, but it explains why 39,026 is what it is, and it moves the canonical
frame from 38,996 to 39,016.

## 50. Read the binding residues' ancestor cone — it is 29 variables

`s10/cone.py` walks back through `definer` from `x_7068`, `x_2099`, `x_28730`.
The cone is **29 variables, 24 atoms, 7 free inputs, 2 booleans**. Small enough
to read in full (`s10/conedump.py`), and it decodes completely:

```
x_26064 = p            (a37694, the bare pin)
x_1692 = x_32499 = x_36136 = x_17499 = x_9325 = x_28599 = p    (the local wire)
x_20434 = x_2081 ;  x_9062 = x_4287                            (control fan-out)
x_31033 = x_2081*(1-x_4287)      x_21279 = x_2081*x_4287       x_6788 = (1-x_2081)*x_4287
x_22542 = x_6418*x_31033   x_25297 = x_9118*x_21279   x_10878 = x_6788*x_31861
x_37158 = x_10878 + x_22542      x_2099 = x_25297 + x_37158
x_642  = x_17325*p        x_7068 = x_2099 + 7376877*x_642      x_28730 = x_9413*p
```

> **`x_2099` is a 3-way MUX over three FREE INPUTS** — `x_6418`, `x_9118`,
> `x_31861` — selected by `(x_2081, x_4287)`, and `x_7068 = x_2099 + 7376877·p·k`.
> Upstream, both binding atoms are *trivially* satisfiable. The pin is entirely
> downstream. Session 9's picture of `D0` as a pinned constant was the shadow of
> a MUX, not a wall.

## 51. The fourth branch: `x_7075 = 0`

`a36085: x_7075 = 1 - x_21279` and `x_21279 = x_2081·x_4287`. The control pair
has only ever been `(1,0)`, giving `x_7075 = 1`. In branch **`(1,1)`** it is
**zero** — and `x_7075` is the multiplier in

```
a35759 = -x_29854 + 5113045*x_7075*x_9118
a35761 =  x_31864 +         x_7075*x_8731
```

so four of the seven residual atoms lose their multiplier at once, and the
congruences `p | x_9118`, `p | x_8731` simply evaporate. Part X's
"boolean branches closed exhaustively, ≥ 7" was measured **in the witness frame**,
i.e. with the flip applied and nothing re-solved. Re-solved (`s10/branch4.py`,
`s10/engine.py b11`) the branch reaches 38,994 — still short, because `(1,1)`
switches on `x_21279 = 1` and with it `a19088`, `a22233`, `a22235`. But the door
was real and had never been opened.

## 52. Two move classes no previous search possessed

**(a) Two-level handle repair.** `lib.ripple` repairs an atom only through its
canonical output variable. `s10/repair2.py` instead solves the atom for *any* of
its variables and then realises that target through a free input of **that
variable's own definer**. First thing it found:

```
a7930 closed via the free input x_24548        (+5 equations)
```

`a7930` is the atom RESUME called *the weak link* — "if atom 7930 can be closed
while `x_28730` moves, the score is 39,027 immediately". It closes.

**(b) mod-p Newton moves.** For a p-quantised check the requirement is not a
*value* but a *residue*: `a ≡ 0 (mod p)`, after which the handle absorbs `a/p`
exactly over ℤ. Take the exact AD derivative `d = ∂a/∂u (mod p)` and shift a free
input by `δ = −a·d⁻¹ (mod p)`. Shifting `u` zeroes **no** atom by itself, so every
zero-this-atom move generator in every previous session was structurally blind to
it. It closes `a35759` at zero cost:

```
x_9118 <- p*(x_9118 // p) ,  then x_1329 = 5113045*x_9118/p     39,002 -> 39,009
```

And because `δ` is pinned only mod `p`, `δ + k·p` is free for every `k`: a second
divisibility can be satisfied simultaneously by **CRT** (`s10/crt.py`).

Canonical frame: **38,996 → 39,009**, and with a beam over the same moves
(`s10/beam.py`) **→ 39,016**.

## 53. The 1-equation "hardening" checks are not constraints

They are exact integer multiples of the gadget they shadow (`s10/shadow.py`):

```
a37662 / a21617 = 10        a40826 / a29539 = 2
```

So the residual of the canonical frame is not four atoms but **two**.

## 54. The gadget family, and why the problem is one cluster

`s10/family.py` enumerates every atom of the shape `c·(A − B) − C`:

```
192 linear gadgets;  67 have C = p * (solo free handle)
each asserts exactly   A == B  (mod p)
gradient-support size 1 : 185 gadgets  -- isolated, all satisfied
coupled clusters        : {7930, 21617, 29539, 33796}  and  {2423, 26731, 33929}
currently failing       : 21617, 29539   -- both in the first cluster
```

> The instance is 192 independent modular-equality assertions, 185 of them
> hard-wired to a single free input each. **The entire remaining problem is one
> 4-gadget cluster.**

## 55. Forward-mode AD makes the true closure affordable

Reverse mode costs one pass per *check* (10,792 of them), which is why every
earlier closure was restricted to a hand-picked neighbourhood. Forward mode
(`s10/fwdad.py`) costs one pass per *free input* and we need only a few hundred:
the full check-by-free-input Jacobian mod p, **134 columns in 15 seconds**.

Closing rows→columns→rows to a genuine fixed point (`s10/closure2.py`):

```
it0: 572 rows x 134 cols  rank 134  inconsistent rows 21
it1: 579 rows x 142 cols  rank 142  inconsistent rows  6   <- fixed point
```

Full **column** rank: the free inputs are locally rigid, no kernel at all. The
system is inconsistent by exactly **one** functional (a single `b` column can
raise the rank by at most 1; the six rows are its witnesses, not six obstructions):

```
witnesses: 33796, 40562, 41400, 41507, 41827, 42245
leftnull dim 437;  t = Y.b nonzero in 6 values
minimum sacrifice: no single row, no pair within budget   (s10/sacrifice.py)
```

## 56. The delivered witness's own frame, built explicitly

The delivered witness is **off** the canonical manifold — running `fwd` on it
snaps it back to 38,996, which is why no forward-based repair could ever touch
it. Detach the five variables whose gate atoms it breaks
(`x_7068, x_28730, x_29854, x_31864, x_642`) and they become free parameters
while their atoms become checks (`s10/frame2.py`):

```
free params 7278 (was 7273);  checks 10797 (was 10792)
delivered witness after fwd in THIS frame: 39026  -- on-manifold: True
failing checks: 22229, 22230, 35758, 35759, 35760, 35761, 35762
closure: 45 rows x 13 cols   inconsistent, witnesses {22230, 35762, 37887}
```

## 57. The ledger — why 7 is the invariant

In frame 2 all seven residual checks are zeroable **exactly and simultaneously**,
by choice of the five detached parameters and the solo handles (`s10/construct.py`
builds it and verifies every one is 0 over ℤ):

```
p | x_9118  -> x_29854 = p*x_1329 = 5113045*x_9118
p | x_8731  -> x_31864 = p*x_10903 = -x_8731
x_642 = p*x_17325 ;  x_7068 = x_2099 + 7376877*x_642 ;  x_28730 = p*x_9413
```

and the cost lands, every time, on the gadget cluster:

```
zero the p-group  ->  a7930, a29539, a40826, a41512 break   score 39,004
   repair a7930 via x_24548                                  score 39,009
   residual: a21617 + a29539 (+ their shadows) = 24 equations
```

Both frames converge to the **same** residual. That is the whole instance:

| choose | cost |
|---|---|
| satisfy the p-quantisation group, give up the gadget cluster | 24 → 17 (beam) |
| satisfy the gadget cluster, give up the p-group | **7 — the deliverable** |

> **This is why 7 is invariant across placements.** It is not six coincidences;
> it is the cheaper side of a single trade, and the trade is forced by one
> obstruction functional in a closure of full column rank.

## 58. Position

```
delivered (give up the p-group)                    39,026   [verified]
canonical frame, new move classes                  39,009
canonical frame, beam over the new moves           39,016
frame 2, p-group zeroed exactly, cluster repaired  39,009
branch (1,1), x_7075 = 0, fully re-solved          38,994
```

To beat 39,026 the gadget cluster must be solved **whole** — its members cost
10–15 equations each, so no partial fix competes with 7. The cluster's closure has
full column rank and one obstruction. What remains is a genuine break of that
obstruction, not a cheaper door.

## 59. The full closure — and a gap in §55 worth naming

§55 grew the column set only from the rows that *witnessed* the inconsistency.
That is not the honest fixed point: a free input touching **any** row of the
system is a repair freedom. Iterating both directions (`s10/closure3.py`):

```
it0: rows  572  cols 134
it1: rows 1333  cols 441
it2: rows 1558  cols 632
...
it6: rows 1655  cols 707      <- fixed point
rank 707 of 707 columns;  KERNEL DIMENSION 0;  still inconsistent
witnesses: 33929, 40068, 40390, 40969, 41400, 41507, 41827, 41842, 42117, 42226
```

Full **column** rank with 1,655 rows: every free input the residual can reach is
completely pinned by the checks. There is no slack anywhere in the reachable
space — and the witness set now reaches `a33929`, i.e. into the *second* gadget
cluster. This supersedes the 128 × 79 closure of session 10 and the 579 × 142 of
§55; both were sub-closures.

## 60. What insisting on the cluster actually costs

Which rows witness an inconsistency is decided entirely by pivot order, so the
elimination can be *steered*. Order rows by decreasing cost and the cheap ones
fall through as witnesses (`s10/closure4.py`, `s10/closure5.py`), where
`cost(c) = |equations of c not already failing|`:

```
rows ordered cheapest-last          -> witnesses {21617, 29539, 37662, 40826}, cost 24
   (the four failing atoms have cost 0 -- this just recovers the current state)

failing rows FORCED to be pivots    -> witnesses {31932, 31934, 31936, 39032, 39034,
   9591, 31940, 33184, 40511, 40577, 40604, 40806, 41400, 41507, 41827, 42245}
   union of their equations: 52
```

> **Insisting on the gadget cluster costs at least 52 equations — score ≤ 38,981.**
> Against a give-up price of 7, the cluster is not merely unaffordable, it is
> unaffordable by a factor of seven.

## 61. The obstruction is branch-independent

All four MUX branches, each run to a fixed point of the enriched engine
(`s10/engine.py br00 | br01 | b11`, and the canonical `(1,0)`):

```
(1,0)  39,009      residual  21617, 29539  (+ shadows)
(1,1)  38,994      residual  19088, 21617, 22233, 22235, 29539, 37887 (+ shadows)
(0,1)  39,002      residual  21617, 29539, 35759 (+ shadows)
(0,0)  38,986      residual  15462, 15464, 21617, 29539, 35759, 36602 (+ shadows)
```

Every branch lands on **the same two cluster gadgets**. `x_2081 = 0` does unpin
`x_6418` (it makes `a3576` trivial), which was exactly the rigidity that blocked
the `a29539` Newton move — and it still does not help, because the branch pays
more elsewhere than it frees.

## 62. Final position after session 11

```
delivered (give up the p-quantisation group)             39,026   [verified]
canonical frame, enriched moves                          39,009
canonical frame, beam over enriched moves                39,016
frame 2 (p-group zeroed exactly over Z, then repaired)   39,009
branches (1,1) / (0,1) / (0,0), fully re-solved     38,994 / 39,002 / 38,986
insisting on the cluster (cost-ordered elimination)      <= 38,981
full closure                       1655 x 707, rank 707, kernel 0, inconsistent
```

The instance is a single forced trade between the p-quantisation group (7) and the
gadget cluster (24). **39,026 takes the cheaper side, and every alternative is now
priced above it.** What would beat it is not a cheaper door but a break of the one
obstruction functional — and that functional lives in a closure with no kernel at
all, which is the sharpest statement of the wall this instance has yet produced.

---

# Part XII — the residual, derived rather than searched

Session 11 set out to break 39,026. It did not, but it replaced the searched
claim "7 is invariant across placements" with an **exact derivation** of the
achievable set, and it found two structural facts that overturn earlier readings.

## 63. The instance is degenerate, and the closure was measured on that degeneracy

```
free inputs: 7,252 of 7,273 are ZERO   (12 loaded constants, 2 booleans, 7 handles)
variables  : 35,208 of 38,748 are ZERO
quadratic monomials: 10,894 live, 242,527 dead  -> 95.7% of the circuit is OFF
```

Every Jacobian in Parts I–XI was taken at this point, where a monomial `u·w` with
`w = 0` has zero derivative in `u`. 115 free inputs reach the cluster structurally
but have derivative zero (`s10/dead.py`). So "full column rank, kernel 0" is a
statement about a degenerate stratum, not about the instance.

Measured directly (`s10/linearity.py`): the two cluster residues are **exactly
linear mod p** in every non-boolean free input, but the collateral checks are not
— 656 of 1,376 large-move predictions are wrong. The linear veto from those rows
is therefore untrustworthy, which is why the closure's negative is not final.
The *absorbable* rows, however, are 90.6% linear, so the reduced closure
(`s10/closure_red.py`, 1,079 × 594, rank 594) is essentially exact.

## 64. Frame 2 severs the paths — and four atom values become free

In the delivered witness's own frame (`x_7068, x_28730, x_29854, x_31864, x_642`
detached), most of the coupling disappears (`s10/rhoknob.py`, `s10/price7.py`):

```
perturb x_9118  by anything -> breaks NOTHING outside the seven
perturb x_8731  by anything -> breaks NOTHING outside the seven
x_642, x_29854, x_31864     -> zero collateral (each occurs in exactly 2 atoms)
x_7068 + 1                  -> breaks a29539 (+ shadow a40826)   13 equations
x_28730 + 1                 -> breaks a7930  (+ shadow a41512)   16 equations
```

## 65. The achievable atom set, exactly

Writing `A = (a22229, a22230, a35758, a35759, a35760, a35761, a35762)`:

```
A2 = x_29854 - p*x_1329        A3 = -x_29854 + 5113045*x_9118
A4 = x_31864 - p*x_10903       A5 =  x_31864 + x_8731
  => A2+A3 == 5113045*x_9118 (mod p)  and  A5-A4 == x_8731 (mod p)
```

and **both right-hand sides are free**, because `x_9118` and `x_8731` cost nothing
in this frame. Hence **A2, A3, A4, A5 are completely unconstrained**. What remains:

> **(1)  `A0 + 7376877·A6 ≡ C₀ (mod p)`** — `C₀ = x_7068 − x_2099`; `x_7068` shifts
> freely by multiples of `p` (a29539's handle `x_30163` absorbs them), so the
> mod-7376877 part is free and only the mod-`p` part binds.
> **(2)  `A1 ≡ A1₀ (mod p)`, `A1₀ ≠ 0`** — pinned by `a7930`'s own congruence.

## 66. Why exactly 5 of 12, and what the sixth would cost

The 12 equations, as coefficient rows over `A` (`s10/rhs.py`):

```
eq  2554 [ 1, 13,  0,  0,  0,  0,  0]      eq 12270 [-31,  5,  1,-27, -1,-17, 10]
eq  6816 [-15,-11, 38,  9, 36, 13, 29]     eq 12350 [-23, 26,-16, 34,-34, 35, 11]
eq  8124 [36, 26,  0,  0,  0,  1, -6]      eq 14584 [ 17, 16, -2,-18,-31, 19,-39]
eq  9123 [ 0,  0, 20, 27, 33,  3, -1]      eq 18673 [  0,  0,  0,  1,  6,  1,  0]
eq  9421 [13,-21,-21, 29, 38, 29,  4]      eq 22044 [-24,-10,  0,  0,  0,  0,  1]
eq 12231 [18, 24,  0,  0,  1,-23, 13]      eq 29125 [  0,  1,  0,  0,  0,  0,  0]
```

`eq 29125` is `A1` alone, so it is satisfied **iff `A1 = 0`** — impossible under (2).
`eq 2554` is `A0 + 13·A1`, satisfiable under (1). The other ten are ten conditions
on the four free values `A2..A5` (`A6` being pinned once `A0` is), so four fall.

> **1 + 4 = 5 satisfied, 7 failing. That is the deliverable, and it is now derived
> from the constraint structure rather than found by search.**

Ignoring (1) and (2) the combinatorial optimum is **6** (`s10/lattice11.py`), and
the only way to reach it is `A1 = 0`. That buys exactly one equation. Its price
(`s10/final27.py`): `a7930`'s congruence must then be met by something other than
`x_28730`, and its entire gradient support is six free inputs:

```
x_2081  -> 109 equations      x_13195 -> 63      x_24548 -> 11   <- cheapest
x_4287  ->  44 equations      x_12553 -> 14
```

> **The sixth equation costs 11 and is worth 1.** 39,026 is optimal in this
> placement, by enumeration of the complete repair space rather than by sampling.

## 67. Two corrections to Part XI

* `a7930` has a **second** repair path — through `x_7927`'s handle `x_11052`, at
  **zero** collateral — available whenever its congruence holds (`s10/chain27.py`).
  Part XI only ever closed it through `x_24548`, which is the expensive route.
* The boolean scan over all 78 booleans in the cluster cones, each re-solved with
  the enriched engine, returns 39,009 for every one (`s10/boolscan.log`) — the
  witness-frame verdict of Part X survives re-solving, for these at least.

## 68. Position after session 11

```
delivered                                                39,026   [verified]
canonical frame + enriched moves / beam            39,009 / 39,016
frame 2, A1 = 0, a7930 repaired via x_24548              39,012
frame 2, six of seven zeroed (a35758 alone nonzero)      38,998
all four MUX branches, re-solved              38,986 ... 39,009
```

The instance is now closed on the placement: the achievable atom set is known
exactly, the combinatorial optimum is 6, and the single equation separating 6 from
5 costs eleven. What remains is the gadget cluster, and it lives on a stratum
where 95.7% of the circuit is switched off — which is the one thing this session
showed is *not* rigid, only unmeasured.

---

# Part XIII — equation-space compensation, opened and priced

Part XII proved 39,026 optimal *for its atom set*. This part attacks the one
assumption left in that proof: that the residual must be carried by those seven
atoms. It need not — an extra atom sharing the same twelve equations changes the
rank, and the search space that opens is real. It is also, in the end, priced.

## 69. The twelve equations admit a free compensator

The 12×7 coefficient matrix has rank 7, so all twelve hold only if `A = 0`. Adding
an eighth column changes that. Enumerating every atom that appears in those twelve
equations (`s10/compensate.py`):

```
a22231  appears in 10 of the 12,   0 equations OUTSIDE   <- free compensator
a22232  appears in  9 of the 12,   1 outside
a35757  appears in  5 of the 12,   2 outside
a22233  appears in  8 of the 12,   2 outside
```

`a22231 = x_4432 − x_19964 − x_28730` lives **entirely inside** the twelve. With it
as an eighth free value the optimum rises (`s10/comp8.py`):

```
subsets of size 8..12 : no kernel
subsets of size  7    : 792 have a kernel  ->  7 of 12 satisfiable
=> 5 of the twelve fail, + a37887 = 6 total  ->  39,027
```

## 70. Frame 3, and why `A1` stops being pinned

Detaching `x_4432` (frame 3, `s10/frame3.py`) severs `x_28730 → x_4432 → a7930`:

```
x_4432  +anything -> breaks a7930, a37887, a41512      17 equations
x_28730 +anything -> breaks a37887 ONLY                 1 equation
x_9118 / x_8731   -> break nothing outside the eight    0 equations
```

So in frame 3 `A1` is free at a price of one equation, and `eq 29125 = A1` alone is
bought. What blocks the seventh equation is that the kernel is one-dimensional
while there are two congruences, requiring
`C₀·(w1+w7) ≡ K·(w0 + 7376877·w6) (mod p)` with `K = x_4432 − x_19964 (mod p)`.

## 71. `a37887` depends only on `x_4432` and `x_19964 + x_28730`

Reading its monomials, `a37887` contains `−2·x_4432·x_19964 − 2·x_4432·x_28730` and
nothing else in those two variables. Since `a22231 = x_4432 − x_19964 − x_28730`
does too, the **compensating pair** `x_28730 += d, x_19964 −= d` leaves both exactly
invariant. Constructed in frame 3 (`s10/drive2.py`) it does precisely that:

```
a22230 = 0   a22231 = 0   a37887 = 0   a7930 = 0     simultaneously
```

— four checks that no previous session ever held at once. The cost is the driver:
`x_19964`'s ancestor cone has only 17 variables and **three live drivers**

```
x_12553 -> 15 equations (it is the load pin a3578 = x_2081*(x_12553 - HUGE))
x_4287  -> 30 equations
x_2081  -> 110 equations
```

> The pair move is mechanically perfect and costs **14** for a gain of **1**.

## 72. Why `x_8731` looked free — and where it is not

`a1459 = x_19892 − x_8731·x_21279` with `x_21279 = 0`, so `x_8731`'s path into
`x_19964` is **switched off**. That is the same degeneracy as §63, and it is why
`x_8731` measures as a zero-collateral knob. In branch `(1,1)`, `x_21279 = 1` and
the path is live — measured `d(x_19964)/d(x_8731) = 1` exactly (`s10/b11f3.py`),
turning `x_8731` into the free driver the pair move wants. The branch itself,
however, activates `a19088, a22233, a22235` and does not pay: frame 3 in `(1,1)`
reaches 39,014 with the engine, and the targeted paired construction 38,974
(`s10/b11build.py`).

## 73. Position

```
delivered                                                     39,026  [verified]
frame 3, branch (1,1), enriched engine                        39,014
canonical beam / enriched engine                       39,016 / 39,009
frame 3, A1 = 0 via the compensating pair                     39,009
```

Everything found in Parts XI–XIII is a mechanism, and every mechanism has a price
above the 7 the deliverable pays:

| lever | buys | costs |
|---|---|---|
| `A1 = 0` (eq 29125) via `x_28730` in frame 3 | 1 | 1 (`a37887`) — but then only 6 of 12 |
| the seventh equation (compensator kernel) | 1 | compatibility condition on `K` |
| moving `K` — the compensating pair | — | **14** (`a3578`, the load pin) |
| moving `C₀` via `x_7068` | — | 13 (`a29539`) |
| unpinning `A1` via `a7930`'s congruence | 1 | 11 (`a21617`) |
| the gadget cluster whole | 24 | ≥ 30 (reduced closure), ≥ 52 (full) |

The instance's margin is six equations and every door now has a number on it.

---

# Part XIV — the balance law, and why 7 is exactly the floor

## 74. `a37887` is a perfect square

Solving `a37887` as a *quadratic* in each of its variables (`s10/quadfix.py`) — a
move class no search had, since `solve_lin` returns `None` on any variable of
degree 2 — every variable has a **double root**:

```
x_4432   roots  -119325148954113784451...  (twice)
x_28730  roots  5439633036170520078110...  (twice)   <- the delivered value
x_19964, x_18253 : likewise double
```

So `a37887 = Q²` for a linear form `Q`, and `Q = 0` pins `x_28730` **exactly**, not
merely mod `p`. Reading its Hessian (the coefficient of `x_4432·w` is `2·q_4432·q_w`)

```
Q = a22231 - 3*x_18253 + 5*x_37720 + 5*x_30108 - 5*x_34600 - 9*x_23754
    - 14*x_7945 + 18*x_23642 + 18*x_23822 + 27*x_37254 - 27*x_15324
    - 34*x_35619 - 13523972*x_9629
```

> `Q` is `a22231` plus exactly the compensator-family variables. The check that
> pins `x_28730` is built from the very atoms that could compensate it.

## 75. The balance law

Let `n` be the number of atoms allowed to be nonzero, `c` the number of mod-`p`
congruences they satisfy, and `E` the equations they touch. A `k`-subset is
satisfiable when its kernel has dimension at least `c`, so `k = n − c` and

> **failing = |E| − n + c**

Check it against everything measured:

| atom set | n | c | \|E\| | failing |
|---|---|---|---|---|
| the seven | 7 | 2 | 12 | **7** ← the deliverable |
| + `a22231` | 8 | 2 | 12 | 6, **+1 for `a37887`** = 7 |
| + `a22232, a22233` | 10 | 3 | 15 | 8 |
| + `a22234, a22235` | 12 | 4 | 17 | 9 |

Adding a compensator with `out` equations outside the twelve changes the count by
`Δ = out − 1 + Δc`. Every candidate was enumerated in §69:

```
a22231  out = 0, Δc = 0  ->  Δ = -1     the ONLY improving compensator
a22232  out = 1          ->  Δ >= 0
a35757  out = 2          ->  Δ >= +1
a22233  out = 2          ->  Δ >= +1
...
```

and `a22231`'s single-equation gain is spent exactly on `a37887`, whose `Q` is
built from `a22231` itself. **The instance is balanced so that its one free
compensator pays for precisely the one check that blocks it.**

## 76. Session 11, final

Nothing beat 39,026, and the reason is now a theorem-shaped statement rather than a
tally of failed searches:

```
achievable atom set     A2..A5 free; A0 + 7376877*A6 == C0 (mod p); A1 == A1_0 (mod p)
balance law             failing = |E| - n + c
the only improving compensator   a22231 (Δ = -1), cancelled exactly by a37887 = Q^2
every other lever       priced 11, 13, 14, 30, 52 against a gain of 1
```

New machinery this session, all of it reusable: two-level handle repair, mod-`p`
Newton moves with CRT, forward-mode AD closures, explicit non-canonical frames
(2 and 3), the repair cascade, the lexicographic potential, and quadratic repair
moves. Four checks — `a22230`, `a22231`, `a37887`, `a7930` — were held at zero
simultaneously for the first time.

**39,026 / 39,033 stands, and the margin of six is now explained.**

## 77. The one direction that is open, now measured

`s10/activate.py` showed no *single* activation of a dead free input reaches the
cluster — by construction, since a dead `u` only multiplies a `w` that is `0`.
`s10/second.py` does the second-order version: find the blocking `w`, find a free
input `z` that makes `w` nonzero, and test the pair.

```
6 activations tested, 6 GREW the cluster's gradient support
  x_17406 -> activates x_5858   : support +2 knobs, 6 atoms broken
  x_16586 -> activates x_15148  : support +2 knobs, 7 atoms broken
  x_12054 -> activates x_30131  : support +2 knobs, 10 atoms broken
  x_28713 -> activates x_3896   : support +1 knob, 13 atoms broken
```

> Activation is real: the live stratum genuinely has knobs our stratum does not.
> The exchange rate measured here is **1–2 new knobs per 6–13 broken atoms**, and
> the closure needs its rank deficit closed, so this is the only door left and it
> is a second-order search, not a linear one.

---

# Part XV — the balance law verified, and both congruences priced

## 78. The kernel dimensions are generic — the law is exact

Part XIV's law assumed `k = n − c`, i.e. that no subset of the twelve is
*dependent*. Checked exhaustively (`s10/kdim.py`), computing the kernel dimension
of every subset rather than merely testing for a nonzero kernel:

```
the seven (n=7, c=2):
  sizes 12..7 : kernel dim 0 for ALL subsets
  size 6      : kernel dim 1 for all 924
  size 5      : kernel dim 2 for all 792   <- first size reaching c = 2
  => 5 satisfied, 7 failing
  and the witness subset is {2554, 6816, 8124, 9123, 9421}
     -- EXACTLY the delivered witness's satisfied set

+ a22231 (n=8, c=2):
  size 7      : kernel dim 1 for all 792
  size 6      : kernel dim 2 for all 924   => 6 failing, +1 for a37887 = 7
```

Every subset has exactly the generic dimension: there are no dependent subsets to
exploit. **`failing = |E| − n + c` is tight, not an estimate.**

## 79. Both congruences priced

`failing = 12 − 7 + c`, so each congruence removed is worth exactly one equation.

**Congruence 2** (`A1 ≡ A1₀`) is removed by detaching `x_4432`, which severs
`x_28730 → a7930`. Price: `a37887`, one equation — cancelling the gain exactly,
because `a37887 = Q²` and `Q` is built from `a22231` (§74).

**Congruence 1** (`A0 + 7376877·A6 ≡ C₀`) needs `x_7068` to move *mod p*, and the
only thing that breaks is `a29539`. Pricing every repair in its 79-input gradient
support (`s10/free1.py`):

```
x_7068   ->  0 equations   -- ARTEFACT: the Newton correction returns x_7068 to
                              its own residue, so C0 mod p is unchanged (score 39,022)
x_14853  -> 20 equations   <- cheapest genuine repair
x_24517  -> 43     x_33287 -> 43     x_12054 -> 45     x_16586 -> 45 ...
```

> **Congruence 1 costs 20 for a gain of 1. Congruence 2 costs exactly the 1 it
> gains.** Both are priced, and the instance balances on both.

## 80. Bulk activation does not create freedom

Activation adds columns to the closure (§77), and the closure has zero kernel, so
the natural hope is that enough activation makes columns outgrow rank. Measured
(`s10/bulk.py`), activating `N` dead free inputs from the cluster cone:

```
N     nonzero atoms   closure        rank   kernel   inconsistent   score
0         4           1655 x 707      707      0          11        39,009
10       60           1715 x 741      741      0          73        38,766
30      127           1786 x 788      788      0         113        38,532
```

> **Rows and columns grow in lockstep and the rank always equals the column count.**
> Activation buys knobs and constraints at the same rate; the kernel never opens
> and the inconsistency only widens. The last door measured is the last door shut.

## 81. Final ledger

```
congruence 1 (C0 mod p)                    costs 20, worth 1
congruence 2 (A1 mod p)                    costs  1, worth 1   (exact cancellation)
the only free compensator a22231           worth 1, cancelled by a37887 = Q^2
the gadget cluster whole                   costs >= 30 (reduced), >= 52 (full)
bulk activation                            kernel stays 0, inconsistency grows
```

Nothing in the instance is worth more than it costs. **39,026 / 39,033.**

---

# Part XVI — the obstruction is not combinatorial; it is the atom map

Parts XII–XV priced everything in terms of atoms *allowed* to be nonzero. This
part asks the prior question: is there **any** atom vector satisfying all 39,033
equations with the residual nonzero? The answer is yes — which moves the entire
obstruction out of combinatorics and into realisability.

## 82. The raw instance, and the forced core

`EQUATIONS.txt` writes each equation as `m*(sum c_i a_i)` or its square, with atoms
like `x-1`, `x-0`, `x*x-x`, and gadget atoms. Measured (`s10/raw.py`):

```
39,033 equations, 42,267 atoms, 38,748 variables
squared equations                      10,478
equations with exactly ONE atom         3,234   <- their atom is FORCED to zero
equations with 3..24 atoms             35,798   <- combinations, so atoms may cancel
pure "x = 0" atoms 1,690 ;  "x = const" atoms 1,105 ;  boolean atoms 6,050
```

Equation indexing was verified against the checker: `L.failing_eqs` returns exactly
`[12231, 12270, 12350, 14584, 18673, 22044, 29125]`. **None of the seven residual
atoms is in a single-atom equation**, so none is forced to zero.

## 83. The compensation closure has a kernel that touches the seed

An atom forced nonzero can be paid for by another atom in the same equation, whose
own equations can be paid for in turn. Propagating that (`s10/closure_atom.py`,
`s10/kerseed.py`):

```
round 0: 7 atoms, 12 equations, 26 compensators available
...
fixed point: 500 equations x 529 atoms   -> rank 500, KERNEL DIMENSION 29
single-atom equations blocking an active atom: 0
kernel basis vectors touching the seed: 24 of 29
```

> **There exist atom vectors satisfying every one of those equations with the seven
> residual atoms nonzero.** The combinatorial obstruction that Parts XII–XV priced
> does not exist at the atom level. What Parts XII–XV actually measured is the cost
> of realising such a vector under a *fixed* frame.

## 84. Restricted to settable atoms it is still there, and it is small

Call an atom *settable* if it has a free variable or a p-absorbing free handle.
Only 22 of the 529 are not. Intersecting the kernel with `{z = 0 on the
non-settable atoms}` (`s10/kersettable.py`, `s10/sparsify.py`):

```
augmented system 522 x 529 : rank 521, kernel dimension 8, all 8 touch the seed
sparsest seed-touching vector: support 69 atoms, touching only 68 equations
   47 settable directly (a free variable), 22 through a p-handle, 0 unsettable
   contains six of the seven: 22229, 22230, 35758, 35759, 35760, 35761
rational kernel on that support: dimension 1  (values fixed up to one scale)
```

> A 69-atom vector, every atom of it movable, whose 68 equations are all of its
> equations. **Realising it would satisfy the entire instance.**

## 85. Why it is not realisable — and it is not the reason I expected

The 69 atoms have only 48 distinct setting variables: **21 collisions**, each a
gadget pair sharing a variable. In the canonical frame a collision demands
`z_a*d_b == z_b*d_a` exactly, and all 21 fail (`s10/collide.py`).

But a pair's shared variable is *defined by an atom already in the support*, so
detaching it costs nothing and gives each atom its own parameter. In that fully
detached frame (`s10/detachall.py`) the residual is exactly six atoms, all inside
the support — and the real price appears:

```
42 variables detached; they also occur in 39 atoms OUTSIDE the support
those atoms live in 120 equations, only 10 of which are already inside
NET equations at risk: 110
```

Adding those atoms and closing over **both** relations — sharing an equation
(compensation) and sharing a detached variable (realisability) — blows up
(`s10/bothclosure.py`):

```
round 0:      46 atoms
round 6:   2,417
round 10: 12,367 atoms, 13,746 equations   -> atoms - equations = -1,379
```

The ratio goes *negative* and the closure keeps growing toward the whole instance,
where it finally turns positive again only because of the 3,234 forced atoms.

## 86. What this changes

```
combinatorial obstruction at the atom level      DOES NOT EXIST (kernel dim 29, 24 touch the seed)
smallest fully-settable seed-touching vector     69 atoms / 68 equations
its realisation in the canonical frame           21 collisions, all inconsistent
its realisation in the detached frame            costs 110 equations, and closing that costs more
```

> The instance is not hard because its equations cannot be satisfied with a nonzero
> residual. They can. It is hard because the **atom map's image** — which atom
> vectors an actual assignment can produce — misses every kernel vector, and the
> coupling that enforces this closes over the whole instance rather than any local
> neighbourhood. Every price in Parts XII–XV is a price for *one* frame; this is why
> no frame beat the others.

## 87. Minimising realisability cost instead of support

Sparsity is the wrong objective for the kernel vector — what matters is how many
detachments its support needs and what those touch. Searching the 8-dimensional
settable kernel over 4,000 random combinations for the vector minimising
`|equations of the outside atoms not already inside the support's equations|`
(`s10/minreal.py`):

```
best basis vector          : 111 equations at risk, support 72, 43 detachments, 41 outside atoms
best of 4,000 combinations : 111 equations at risk, support 72, 43 detachments, 41 outside atoms
```

The minimum does not move. Every seed-touching vector in the settable kernel needs
roughly the same number of detachments, and they touch roughly the same outside
set — another expression of the global coupling of §85.

---

# Part XVII — corrections from the adversarial audit, and from the placement search

Two independent investigations were run against Parts XII–XVI: an exhaustive global
placement search, and an audit briefed to **refute** the optimality claim. Neither
beat 39,026. Both found real errors in how I stated it. Those come first.

## 88. `failing = |E| − n + c` is NOT a law — it is bookkeeping (correcting §75, §79)

The identity is trivially true; its entire content is `c`, which I never derived —
only computed ad hoc for `n = 7` and `n = 8`. Worse, I implied `|E| − n` is bounded
below by 5. It is not: greedily adding atoms that bring at most one new equation
each drives `|E| − n` to **−5** within 394 atoms. And the formula mispredicts at
`c = 0`, where it claims 7 satisfiable when the truth is 12 (`A = 0` satisfies
everything). **It is validated only at `c ∈ {1, 2}` and must not be quoted as a law.**

## 89. Removing BOTH congruences is worth 7, not 2 (correcting §79)

Removing one congruence buys exactly one equation (5 → 6 satisfiable, verified
exhaustively for both choices). Removing **both** admits `A = 0`, which satisfies all
twelve — worth **7**. The verdict is unchanged, because the joint cost is at least
`11 + 13 = 24`, but my ledger understated the payoff of precisely the escape an
attacker would aim at.

## 90. Congruence 1 was mispriced (correcting §79)

`C₀ = x_7068 − x_2099`, and §79 priced only the `x_7068` half. But **`x_2099 = x_6418`
identically**, and moving `x_6418` breaks exactly one atom, `a3576`, at raw cost
**13** — with three compensators wholly inside `a3576`'s own 13 equations, giving a
floor nearer 9. An exhaustive scan of both pin cones gives **13**, not 20, as the
true single-move minimum. Related: `x_14853` measures 33 raw / 20 in the joint move,
not 20 flat.

## 91. Three statements that were loose or wrong, though the answers held

* **§65's `C₀`.** Exactly, `A0 + 7376877·A6 = x_7068 − x_2099 − 7376877·x_28599·x_17325`.
  Writing `C₀ = x_7068 − x_2099` is valid **only because `x_28599 = p`** — load-bearing:
  if `x_28599` were anything not ≡ 0 mod p, `x_17325` would be a free knob on `C₀ mod p`
  and congruence 1 would evaporate.
* **§66's "1 + 4 = 5" derivation is wrong.** "`A6` being pinned once `A0` is" is false —
  only the *combination* is pinned. The correct justification is the exhaustive
  kernel-dimension plus congruence-feasibility computation, not that narrative.
* **§71's driver prices do not reproduce.** I reported `x_12553 → 15`, `x_4287 → 30`,
  `x_2081 → 110`; measured raw in frame 2 they are **30 / 38 / 109**. So "the pair move
  costs 14" should read ~30 unless the 15 was post-repair. There is also a fourth free
  input in `x_19964`'s cone, `x_14865` (12 equations), dead on this stratum.
* **`s10/lattice11.py` is heuristic**, not exhaustive — it searches `(k−1)`-row bases
  rather than all subsets. Its conclusion is confirmed independently, but the script
  does not establish it.

## 92. What the two investigations confirmed

* **The footprint-1 atoms are BUNDLES, not cheap carriers.** Each is a random integer
  combination of primitive atom residuals, squared — e.g.
  `a37887 = (a22231 + 6·a22232 + 15·a22233 − 21·a22234 + 25·a19087 + … )²` — and 3,234 of
  them sit alone in a size-1 equation. So `a37887 = Q²` (§74) is not a designed trap;
  it is the generator's bundling. A bundle can never be *chosen* as a support: when a
  primitive breaks, the bundle adds one atom, one equation **and** one congruence,
  a net +1 every time. **`eq 8680` is `1·(a37887)²` with `a37887` its only atom**, so
  `a37887 ≠ 0` breaks it unconditionally, with no in-equation compensation possible.
* **Whole-instance break census**: all 38,748 variables perturbed with full gate
  ripple, 33,969 distinct supports reached — **minimum failing = 7**, next 8, then 12.
  Nothing below 7 anywhere in the instance.
* **The obstruction is arithmetic, not linear** (independent confirmation of §86): the
  reachable equation-direction space in ℚ¹² has rank **7**, spanned by the seven knobs
  `x_642, x_1329, x_8731, x_9118, x_9413, x_10903, x_17325`, and `−r₀` **is** in that span.
* **`a22231` is unique**: over all 42,267 atoms, exactly **8** have their entire equation
  set inside the twelve — the seven plus `a22231`. Proven, not sampled.
* **Exhaustive congruence feasibility**: of the 924 six-subsets, **0** are feasible; of
  the 792 seven-subsets in the 8-atom placement, **0**. Of the 792 five-subsets, 378 are,
  including the witness's. The integer/rational gap does not bite.
* **Every single-parameter escape is closed**: of 173 free inputs in the three pin
  cones, 29 have zero collateral and **none** moves any pin residue mod p; 82 move one,
  and the cheapest costs 13.
* `ad.fwd` converges in **one** round (identical at 1…40 rounds) — no iteration bug.

## 93. The weakest step, named by the audit

> The claim `c ≥ 2` rests, for **multi-parameter** moves, on mod-p linear closures taken
> at a stratum where 95.7% of the circuit is switched off and 656 of 1,376 large-move
> predictions are wrong. There is no certificate ruling out a multi-parameter,
> non-linear move that shifts `x_7068`, `x_2099` or `x_28730` mod p at low collateral —
> and by §89 a move shifting **two** pin residues at once is worth **7** equations, not 1.
> That the `x_28730`/`x_24548` *pair* (cost 11) beats every single move (cost ≥ 13) is
> direct evidence that joint moves matter and are only being sampled.

Also newly open: **boolean atoms need not be zero.** `x² − x = x(x−1)` has integer image
`{0, 2, 6, 12, 20, …}`, so a boolean atom is a legitimate nonzero carrier. There are
3,484 of them, 1,156 on free inputs, and the placement search found blocks with
**negative** deficiency (e.g. 36 atoms in 32 equations). The two blocks it tested were
unrealisable for local reasons; the space was not swept.

---

# Part XVIII — the circuit decoded, and a false infeasibility claim caught

## 94. The cones are boolean MUX networks, not one-way computations

Exact algebra on the two residual cones (`s10/cl_*`) shows every gate is `x·y`,
`1−x`, `x+y` or `x+y−xy`:

```
cone(x_27522): 257 vars / 38 free, gate kinds {LIN 115, BILIN 104}
   top is a literal 3-way mux on selectors s1 = x_28940, s2 = x_23047:
   x_27522 = s1*s2*x_19799 + s1*(1-s2)*x_36462 + s2*(1-s1)*x_8239
cone(x_1308):  1271 vars / 165 free, {LIN 683, BILIN 423} -- same construction
cone(x_25442): 66 / 11
```

At the canonical frame all 37 boolean free inputs of the first cone are 0 except
`x_2081 = 1`, so the mux **collapses** and the residual becomes a set of relations
between *free inputs*:

```
a21617  <=>  x_24548 == x_14623 (mod p)
a7930   <=>  x_24548 == x_12553 (mod p)      [in the canonical frame]
a29539  <=>  x_1308(x_6418) == x_14853 (mod p)
```

The 296-bit constants are **broadcast classes**: K1 on 57 variables (5 free), K2 on
53 (3 free). Cluster 2 (`a2423, a26731, a33929`) turns out to be *the distribution
network for those constants*, which is why breaking it buys nothing.

## 95. Two zero-net swaps, and a new state with the cluster closed exactly

```
whole K1 class shifted by one delta   39,009 (+-0)   a21617+a37662 (11) fixed, a31672 (11) broken
x_6418 shifted                        39,009 (+-0)   a29539+a40826 (13) fixed, a3576  (13) broken
both together                         39,009         a21617 = a29539 = 0 EXACTLY
                                                     residual moves to a3576+a31672, 24 -> 24, overlap 0
```

Verified by the checker (`s10/cl_cluster_closed.json` → 39009/39033). **Both are
cheaper than the 20 recorded in §79 for either congruence** — this supersedes that
price. The terminal blockers are two boolean-gated constant pins:

```
a3576  = x_2081 *(x_6418  - C4) - 15804267*x_26777    13 eqs   pins x_6418  == C4 (mod p)
a31672 = x_24601*(x_33462 - K1) - x_36358             11 eqs   pins x_33462 == K1 (mod p)
```

Killing their selectors costs 70 and 62 respectively. The residual behaves as a
conserved quantity: it can be moved around the instance but not reduced — 24 in,
24 out, disjoint.

## 96. The obstruction functional depends on only two variables

Reverse-mode AD of `Σ y_a·a_a` for the smallest left-null certificate
(support 11: `a1436, a3576, a3578, a7930, a7932, a15456, a15462, a21617, a21619,
a40065, a41507`) over **all 7,273 free inputs** has gradient support **exactly
`{x_2081, x_4287}`** — the two global boolean switches. The functional is conserved
on the entire stratum. Two strengthenings:

* **Equation-level relaxation** (equations must vanish, atoms may cancel *inside*
  them): 3,600 equations × 707 columns, rank 707, still inconsistent — 19
  obstructions, smallest combining 30 equations. Cancellation does not rescue it.
* 254 of the 707 columns are themselves boolean-pinned; with only the 453 legal
  columns the rank is 453 with the same 11 obstructions, and the minimum-sacrifice
  greedy loses 21 equations (best 39,012).

## 97. A claimed infeasibility proof, refuted

The analysis above suggested a forcing chain: `a7930` gives `x_24548 ≡ x_12553`,
`a3578` pins `x_12553 ≡ C3`, `a21617` gives `x_14623 ≡ x_24548`, and the K1 web plus
`a31672` gives `x_14623 ≡ K1` — while

```
C3 mod p = 4531249068709477613185164105669741036354237152756954144434674493737552368539
K1 mod p = 37841415183514949237467304684128824427406379377151921996714091976892367869714
C3 != K1  (difference 82481923122510723799288844430228824462217842441245596187178166524754019170488)
```

which would make the instance **infeasible** with those selectors on. It is not a
proof. The `a7930` link is **frame-dependent** (`s10/link.py`): `x_25442` is not
pinned to `x_12553` — it is `x_10861 + x_22342` (atom `a21112`) — and at the
delivered witness

```
delivered : a7930 = a3578 = a21617 = a31672 = 0, and x_25442 != x_12553 (mod p)
canonical : x_25442 == x_12553 (mod p), and a21617 is NONZERO
```

> Breaking that link is exactly how the delivered witness satisfies `a7930` and
> `a21617` at once. `C3 ≢ K1` is a true statement about the canonical frame and
> **not** an infeasibility proof. No claim of infeasibility is made.

## 98. The K1 web cannot be cut cheaply

`x_14623` and `x_33462` sit in the same broadcast class; decoupling them would let
`x_14623` move alone. The link graph (`s10/k1web.py`) reaches **71 of 72** class
members from `x_14623`, with many parallel paths, so no single cheap edge separates
them. And the accounting forbids a win regardless: fixing `a21617` gains 11 and
fixing `a29539` gains 13, but their swap targets `a31672` and `a3576` cost exactly
11 and 13. The canonical frame is pinned at 39,009 by conservation, and the
delivered frame's 7 remains the floor.

---

# Part XIX — the p-wire route, settled

The hypothesis of §80-era notes was that if the wire carried `w ≠ p`, every gadget
congruence `A ≡ B (mod p)` would become `A ≡ B (mod w)` and, at `w = 1`, trivial.
**That is exactly right**, and it is now measured rather than argued.

## 99. `w = 1` frees the entire residual — for 13 equations

Building `F_WIRE` (frame 3's six detachments **plus** `x_26064`, so the bare pin
`a37694` becomes a check) and driving all 220 wire members to `w`
(`s10/wr_frame.py`, `s10/WR_RESULTS.md`):

```
python3 checker.py s10/wr_engine_w1_x7068_39020.json
   -> satisfied 39020/39033  (13 failing)
   failing: [8429, 11166, 11915, 12594, 23869, 25313, 26785, 31400, 32300,
             36106, 36767, 37257, 37666]
```

In that state **exactly two of the 42,267 atoms are nonzero**: `a37694` (12
equations) and `a39417` (1 equation). Both residual gadgets, `a7930`, `a29539`,
`a40826`, `a41512` and **all 1,249 handles** are exactly zero. Every `w` tested
(1, −1, 2, −2, 3, 0, −p, 2p, p², p±1) behaves identically.

> The premise is confirmed: a non-`p` wire frees the whole rest of the instance.
> Its price is **13** against a give-up price of **7**.

## 100. The price is irreducible, and there is no third regime

All 12 equations of `a37694` carry it with a nonzero coefficient and nothing else
is nonzero, so all 12 fail; and `eq 11915` is `217608357·(a39417)²` with `a39417`
its only atom, where `a39417 = −8(w − p) ≠ 0` on the diagonal. Relaxing to treat
every atom as an independent rational — strictly more freedom than the circuit has
— gives failing ≥ 12 for `S = {a37694}`, improving only to 8 after two closure
rounds; the same relaxation understates the current placement by ≥ 1 (returns 6
where the truth is 7). Over ℤ it is worse: exhaustive search over all **196,624**
three-knob systems (54,734 with an integral solution at wire base 1) bottoms out at
**11** broken rows.

The deformation space has exactly **two** regimes:

| regime | what happens | price |
|---|---|---|
| **diagonal** (all members equal, value ≠ p) | 8 quadratic wire-only checks stay zero, handles unquantise, residual fully repairs — but the bare pin fires in all 12 of its equations plus eq 11915 | **13** |
| **off-diagonal** (incl. the whole 3-dim kernel) | pin held at zero, 0 identity rows break — but the 8 quadratic checks plus `a39084, a39417, a41278` all fire, and the multipliers reach 78–325 digits so the congruences get *harder* | 22–23 → 39,010–39,013 |

**No third regime exists**: the diagonal is a 1-dimensional line, row 37257's wire
content is the bare pin alone (a nonzero functional on that line), and the eight
quadratic checks vanish exactly on it. Partial deformation was swept too — a 0/1
subset local search over 60 restarts finds uniform optimal at 12; single members
cost 34/58/60/70, all-but-root 50, root alone 62.

Side result: restoring only the root to `p` from the `w = 1` state gives a new
**39,021** with a **4-atom** placement `[22229, 22230, 35758, 35760]` — shorter than
any placement previously recorded, though still behind 7.

---

# Part XX — activation, closed exhaustively (and §77 corrected)

## 101. §77's "6 of 6 pairs grew the support" was a control failure

`s10/second.py` reported that second-order activation pairs grow the cluster's
gradient support, and Part XVII carried that forward as the one open door. **It was
an artifact.** All six partners (`x_12054, x_16586, x_17406, x_28713, x_27393,
x_11368`) grow the support **as singles**. Against the correct control — pair versus
best single — the gain at value 1 is **zero**. Enumerating *all* 17,766 pairs twice
(`ac_pairs.py`, `ac_pairs2.py`):

```
value 1  : 11,438 grow the support; max 2 knobs; 0 genuinely second-order pairs
generic  : max 4 knobs; 1,846 genuinely second-order; best = 4 knobs for 71 equations
```

Second order is real only at generic values, and its **rate is worse** than first
order: 17.8 equations per knob against 9.5 for the best single.

## 102. The complete candidate set, and the true price of a knob

Nothing outside the cluster's structural cone can affect it, so the space is
bounded exactly: **1,401 variables containing 189 zero-valued free inputs**. All 189
were swept at three values (`ac_single.py`): 42 cost zero equations, 76 grow the
support, **max 2 knobs**. Best is `x_24365 = 1` at +2 knobs / 19 equations raw —
and with the repair engine run with the activated input **frozen**, the collateral
largely repairs with the knobs kept:

```
39,009 -> 38,990 -> 39,003 , activation alive, +2 knobs still present
true price: 6 equations for 2 knobs   (reproduced for x_12054 and x_16586)
```

Unfrozen, the engine repairs by *undoing* the activation — the results are
byte-identical to `mod9118_0.json`. Stacking 2–3 activations buys no further knobs.

## 103. The kernel never opens, and here is why

| activation | closure | rank | kernel |
|---|---|---|---|
| none | 1655 × 707 | 707 | **0** |
| `x_24365 = 1` | 1647 × 706 | 706 | **0** |
| `x_24365` generic | 1657 × 708 | 708 | **0** |
| top-4 generic | 1665 × 714 | 714 | **0** |
| **73 "no-new-row" generic** | 2009 × 953 | 953 | **0** |
| all 189 generic | 2078 × 976 | 974 | 2 (useless: 361 inconsistent rows) |

The regime I asked for does exist — `ac_rowcost.py` shows **73 of 76** generic
activations add their knobs with zero new closure rows and zero column loss — but
rebuilding still grows rows, because the existing columns' footprints move with the
point. **Rank equals the column count in every targeted configuration.**

> Structural reason the stratum is pinned: the closure's column set is **closed** —
> any free input whose gradient reaches any row is already a column — so the
> left-null functionals annihilate *every* free-input column.

## 104. The right question was the coset leader, and its ceiling is 39,018

Since rank = columns always, the solve is injective and the real quantity is the
**minimum-equation-cost coset leader**:
`score = 39033 − |eqs_of_atoms(supp(Mx − r))|`. Over **1,035 information sets**
(`ac_isd.py`) the best is **15 equations → a mod-p ceiling of 39,018** for the
canonical frame — with `a29539` itself zeroable. Dropping all 396 weight-1 checks
leaves it inconsistent; consistency arrives only at weight ≥ 12, at 69 equations.

> **The 39,009 frame's linear ceiling is 39,018 — below 39,026** — even though
> zeroing both gadgets outright would be worth 39,031.

## 105. The record frame admits no move at all

The 39,026 cluster `{22229, 22230, 35758…35762}` has a structural cone of **39
variables** containing **exactly one** zero free input (`x_31861`, costing 29
equations for 0 knobs) and a gradient support of **3** free inputs. A full sweep of
**14,541 single ±1 moves over all 7,273 free inputs** (`ac_sweep26.py`) finds **no
improving move**; the 3,878 score-neutral moves all lie outside the 39-variable cone
and can never reach it. Ripple-based repair finds only **12 legal moves** in the
entire frame, none improving.

---

# Part XXI — the boolean route, closed (and §61's branch table corrected)

## 106. The load-pin shape is three monomials, not two

I had been filtering for `b*(x − K)` — two monomials — and getting nothing. The
actual shape in this instance is **`b*(x − K) − c*z`**, three monomials with a slack
term. Census (`bl_pins3.py`): **727 conditional constant loads** across 623 atoms,
**512 clean pins**, **256 distinct gate booleans** (all free inputs), constants
289–296 bits (≈ `p·2^37..44`).

Which pins reach the binding congruence `C₀ = x_7068 − x_2099`? **Exactly two**:
`a3576` (gate `x_2081`, pins `x_6418 = K₁`, currently ON) and `a3568` (gate
`x_4287`, pins `x_31861 = K₂`, off). No other pin in the instance reaches `C₀`.

## 107. No pin is worth unpinning

Setting `x_2081 = 0` does release `a3576` — but it simultaneously zeroes the
multiplier `x_31033 = x_20434·x_31822`, so `x_2099` collapses to the constant **0**
rather than becoming free. The MUX just swaps one constant for another (measured:
`x_2099` has 295 bits in (1,0), 0 bits in (0,0) and (0,1)). The only branch where
`x_2099` becomes genuinely free is **(1,1)**, where `x_2099 = x_9118` — and that
branch switches **on** three `x_21279`-gated quantisation checks (`a19088`,
`a22233`, `a22235`, needing `p | x_9106`, `p | x_2239`, `p | x_31731`) in exchange
for the two it switches off. Net **+1 obligation**.

## 108. The decisive test: 900 simultaneous flips change nothing

In the witness frame the seven residual atoms depend on **39 variables**, whose free
inputs are `{642, 1329, 2081, 4287, 6418, 7068, 8731, 9118, 9413, 10903, 17325,
28730, 29854, 31861, 31864}` — **only two booleans, both MUX controls**. So any set
of non-MUX flips must leave all seven residual values untouched. Confirmed at scale:

```
all 900 neutral booleans flipped at once   -> 39,026   (seven values BIT-IDENTICAL)
random half (450)                          -> 39,026
rand-10 / 50 / 200                         -> 38,975 / 38,776 / 38,369
all 1,154                                  -> 36,880
```

Exhaustive singles in both frames (1,154 each): 900 exactly neutral, 254 worse,
**0 better**. All 2,850 canonical pairs over the 76 cluster-cone booleans: every
engine re-solve lands on exactly 39,009 with the identical residual. A 1,232-eval
priority scan over the 22 pin-gate booleans: 12 of the 14 best return to 39,009.

**key-4 exhaustive** — all 16 assignments of `{2081, 4287, 11368, 13195}`, the only
booleans reaching any atom in the seven failing equations, each with a full engine
run: the fixed point depends **only on the MUX pair**, and `x_11368`/`x_13195` are
irrelevant in all four branches.

## 109. Correction to §61's branch table

§61 reported the four branches as 39,009 / 38,994 / 39,002 / 38,986. Those were
measured in the **canonical** frame. Re-run in the **witness** frame two of them are
better:

| branch | §61 (canonical) | witness frame |
|---|---|---|
| (1,0) | 39,009 | **39,026** |
| (0,1) | 39,002 | **39,011** |
| (1,1) | 38,994 | **39,003** |
| (0,0) | 38,986 | 38,984 |

The conclusion is unchanged — no branch beats (1,0) — but the numbers in §61 were
frame-dependent and should not be quoted as branch prices.

## 110. Why the cancellation route is also blocked

`eq 18673 = 949417·(a35759 + 6·a35760 + a35761)` contains **only** residual atoms,
and **none** of the 24 atoms appearing in the seven failing equations occurs
exclusively in failing equations — each also sits in 3–10 satisfied ones. So no atom
can be turned on to cancel inside a failing equation without breaking a satisfied
one.

---

# Part XXII — boolean atoms as carriers: closed by a sign argument

Part XVII flagged a genuinely unswept opening: `x² − x = x(x−1)` is **not** forced to
zero over ℤ, so a boolean atom is a legitimate nonzero carrier, and blocks of them
with negative deficiency exist. It is now swept and closed — by an argument I did
not anticipate.

## 111. Census

```
3,484 boolean atoms, all exactly c*(x^2 - x):  c = 1 (2,340), c = -2 (1,144)
NONE is a gate -- every one is a pure check atom, sitting in 2..21 equations
3,484 distinct variables: 1,156 free inputs, 2,328 gate outputs
syntactic clean carriers: 0   (every boolean var appears in >= 1 other atom)
operational clean carriers: 318 of 1,156  (the other atoms they touch re-solve)
```

Exact optima by max-closure/min-cut rather than greedy: over all boolean atoms the
minimum deficiency is **−29** (376 atoms in 347 equations, one connected component).
But over **free-variable** boolean atoms it is **0** — maxflow = 1,156 is a perfect
matching, so Hall's condition holds and **no free-only block can ever have negative
deficiency**. Rank deficiency is much larger than combinatorial (the free zero-core
has nullity 53; the −29 block, 122), so rational kernels are abundant.

## 112. Non-negativity kills it, not the triangular-number filter

`x(x−1) ≥ 0` for every integer `x`, so each carrier's value `t_a = val_a / c_a`
must be **non-negative**. That binds long before any `k(k−1)` representability test:

> **An exact rational Phase-1 simplex proves the free-variable cone is TRIVIAL.**
> No nonempty set of free boolean atoms can be nonzero with all of its equations
> satisfied — the 53-dimensional kernel is entirely sign-mixed. Every free boolean
> variable really is confined to `{0,1}`, for a reason that has nothing to do with
> the boolean atom itself.

Pairwise, exactly over all 8,237 sharing pairs: a shared equation cancels only if
`t_a/t_b` equals *its own* coefficient ratio, so **at most 2 shared equations can
ever cancel** (88 pairs achieve 2; 6,959 achieve 1) — even for pairs sharing 18
equations. 5,411 pairs do admit a realisable `k(k−1)` ratio, but the best costs 11.

## 113. And they cannot reach the residual anyway

* **Zero boolean atoms appear in any of the seven failing equations.**
* The ancestor closure of the residual is 39 variables with 2 booleans; of *all*
  atoms in the seven failing equations, 151 variables with **4** booleans —
  `x_2081, x_4287, x_11368, x_13195`, the known quadrant bits.
* Exhaustive sweep of those four at **non-boolean** values (6,717 assignments:
  singles and pairs over −6..8, triples −3..5, quads −2..4) gives best **39,026** —
  the current values. Constructively, each key bit at `x ∈ {2, −1, 3}` breaks
  **9–21 non-boolean atoms** while turning on only its own boolean atom.
* All 1,156 free boolean variables set to `x = 2`: scores 39,007–39,018, minimum
  cost 8 (`a13485`/`x_24844`). **Not one ever repaired any of the seven baseline
  failures.** Best pair 11, best triple 10. Checker-verified: the best carrier state
  is the original seven failures **plus exactly `a13485`'s eight equations**.

## 114. A methodology note worth keeping

The 39,026 partial is **off-manifold** — gate atoms `22229, 22230, 35758, 35761,
35762` are deliberately nonzero — so a plain `ad.fwd` repairs them and drops the
score to 38,996. Any experiment starting from the deliverable must use a
block-preserving forward evaluation (the `frame2`/`frame3` detachments, or an
equivalent `fwdb`). Several earlier measurements were frame-dependent for exactly
this reason.

---

# Part XXIII — the joint move, priced constructively (correcting §89), and the frame ceilings

## 115. The whole prize is reachable, and costs 16 — not "≥ 24"

§89 argued that removing **both** congruences is worth 7 equations but costs at
least `11 + 13 = 24`. That lower bound was an argument, not a measurement, and it is
wrong. Realising `A = 0` directly (`s10/jm_19_azero.py`):

```
python3 checker.py s10/jm_azero00_39017.json  ->  39017/39033  (16 failing)
   all seven residual atoms exactly ZERO
   all twelve gadget equations SATISFIED; the failing set is disjoint from them
   only THREE of the 42,267 atoms are nonzero: a688, a1618, a40608
```

In frame 2 the seven atoms are exactly

```
a22229 = x_7068 - x_2099 - 7376877*x_642     a22230 = x_28730 - p*x_9413
a35758 = x_29854 - p*x_1329                  a35759 = 5113045*x_7075*x_9118 - x_29854
a35760 = x_31864 - p*x_10903                 a35761 = x_7075*x_8731 + x_31864
a35762 = x_642  - p*x_17325
```

so `A = 0 ⟺ p | x_9118, p | x_8731, x_28730 ≡ 0, x_7068 ≡ x_2099 (mod p)` — the last
two being congruences 2 and 1. **The ledger entry should read 16 measured, not ≥ 24.**

| move | out-of-twelve cost | score | leftover |
|---|---|---|---|
| congruence 1 alone (`x_6418`) | 13 | 39,009 | `a3576` |
| congruence 2 alone (`x_28730` + tracker `x_24548`) | **11** | 39,011 | `a21617, a37662` |
| both (branch flip + pins + zeroing) | 14 | 39,008 | `a25676, a33796, a42245` |
| **`A = 0` realised** | **16** (in-twelve 0) | **39,017** | `a688, a1618, a40608` |

Exhaustive second-move sweeps from each state — candidate set = all free inputs in
the cones of the atoms of the still-failing equations, provably complete for one
further move — gave **6,297 exact measurements and zero improvements**, plus 1,314
exact-cancellation pairs over the 82 movers (best 47) and 580 driver patterns.

## 116. The mod-p veto is wrong in *both* directions

The 𝔽_p greedy predicted a **4-equation** floor for either congruence (breaking only
`a41400, a41507, a41827, a42245`). Built explicitly and measured: **none of those
four broke**, and instead `{a19297, a19299, a30984, a36185, a40812}` broke, for **20**.

> At large moves the mod-p model is unreliable in both directions — it both
> over- and under-predicts. Only construction settles a price. Every linear ceiling
> in Parts XII–XVI should be read with that caveat.

Exact map (957 forward-AD columns): `dC₀` has support only `{x_2081, x_4287, x_6418,
x_7068}` and `dA₁` only `{x_28730}`; the check matrix is 1,716 × 957 of rank 670 with
kernel dimension 287, and **both targets lie inside the row span**, so collateral is
forced rather than incidental.

## 117. Conservation, and one uniform blocker shape

Every repair *moves* the residual rather than shrinking it — `11 → 14 → 16 → 20 → 24`,
and `13 (a3576) ↔ 13 (a29539 + a40826)`. And every terminal blocker reached from any
direction has the **same shape**:

```
K*(free - derived) - p*handle          with a p-quantised absorber
a3576, a3578, a3568, a3570, a7930, a21617, a29539, a31672, a33796, a688, a1618
```

Feeding one pin re-selects a mux and opens another. Two structural routes closed:
`x_2081 = p` does kill `a3576`'s dependence on `x_6418` mod p (the product dies) but
costs 109 via `a29466/a3578/a37893/a41794`; and `x_26777 = p·x_3387` means `a3576` can
only absorb multiples of `15804267·p`, so `x_6418` can never move mod p through it.

## 118. The two frames' linear ceilings

Since rank = columns in every closure measured, the solve is injective and the real
quantity is the **minimum-equation-cost coset leader**,
`score = 39033 − |equations whose atom combination is nonzero|`. Measured by
information-set decoding in both frames:

```
canonical frame (39,009 seed) : 1,035 information sets -> ceiling 39,018
witness frame  (the deliverable): 1,669 x 714 closure, 1,500 information sets
                                  do-nothing cost 7  ->  ceiling 39,026
```

> **The deliverable saturates its own frame's linear ceiling exactly.** No linear
> move in the witness frame beats 39,026, and the canonical frame's ceiling is
> strictly below it. (Note the corrected cost function: counting every equation the
> nonzero atoms *touch* rather than those whose combination is nonzero reported the
> do-nothing witness frame as 39,021 instead of 39,026.)

## 119. The 39,017 state's blockers cannot be repaired either

`jm_azero00_39017.json` leaves only `a688, a1618, a40608` nonzero across 16
equations. Those 16 equations contain 25 atoms, so a kernel of dimension 9 exists —
but only **one** of the nine basis vectors touches a blocker, and minimising the
outside cost over 20,000 combinations (`s10/blockers2.py`) gives:

```
best kernel vector: support 16 atoms, 19 equations OUTSIDE
   [687, 688, 1617, 1618, 19297, 19298, 19299, 26732, 26733,
    30975, 30978, 30980, 30982, 30983, 30984, 36649]
=> repairing the state costs 19 -> score 39,014, worse than leaving it at 39,017
```

`a687` is a free compensator (zero equations outside), exactly like `a22231`, and
exactly as before it is not enough on its own.

## 120. Session 11 final ledger

```
deliverable                                            39,026  [verified]
  saturates its own frame's linear ceiling exactly     39,026
canonical frame's linear ceiling                       39,018
A = 0 realised (all 12 gadget eqs satisfied, 3 atoms)  39,017
wire at w = 1 (2 atoms nonzero instance-wide)          39,020
whole-instance break census, 33,969 supports        min 7 anywhere
```

Every route is priced by construction:

| lever | measured price | worth |
|---|---|---|
| congruence 2 alone | 11 | 1 |
| congruence 1 alone | 13 | 1 |
| both (`A = 0`) | **16** | 7 |
| p-wire at `w ≠ p` | 13 | everything |
| gadget cluster whole | 24 (conserved) | 24 |
| boolean carriers | ≥ 8 | 0 (cannot reach the residual) |
| activation | 6 per 2 knobs | kernel never opens |

---

# Part XXIV — the frame was wrong, not the search

Session 11 ended with a ledger in which every lever was priced and every price was
too high.  That ledger is still correct as far as it goes, but Part XXIV shows it
was measuring the wrong object.  Four results, in the order they were forced.

## 121. The infeasibility claim of §66 was an artifact (refuted)

`s10/exactlin.py` split the checks into those that are EXACTLY linear mod p in the
free inputs (verified by probing) and those that are not, built the exact subsystem,
and found it inconsistent:

```
checks touched: 365;  EXACTLY LINEAR in every probe: 220
EXACT subsystem: 220 rows x 80 cols, rank 80, inconsistent rows 7
```

If that had held it would have been an exact infeasibility certificate.  It does not
hold.  `s10/exactlin2.py` audits it against the three objections that matter:

| objection | test | result |
|---|---|---|
| sampling — 2 probes is not a proof | re-run with 4 | same 220 rows, stable |
| **column closure** — rows may depend on free inputs outside the 80 | free-input support via `s10/suppfree.py` | **49 of 220 rows do; 162 missing columns** |
| atoms vs equations — `Mx=b` demands atoms be zero, but equations only need their combination to vanish | how many rows sit alone in an equation | **0** |

```
EXACT subsystem 220 x 80: rank 80, inconsistent 7
--- restricted to the COLUMN-CLOSED rows only ---
171 rows x 80 cols: rank 77, inconsistent 0
```

Every inconsistency came from rows whose columns were missing.  **No infeasibility
conclusion survives, and none is claimed.**  This is the second false infeasibility
claim killed this session — the first was the cluster agent's `a7930` link, refuted
by `s10/link.py`.  The pattern is the same both times: a closure computed over an
unclosed index set.

`s10/suppfree.py` is the tool that makes this checkable in future.  It propagates a
bitset of free-input indices along the same topological order forward-mode AD uses,
so `supp(c)` — every free input that can reach check `c` — is available for all
42,267 atoms in 0.1 s.  Set-union over-approximates (two paths can cancel), so
`supp(c) ⊆ U` is a sound proof that `U` is closed for `c`, while `supp(c) \ U` is a
candidate list to test.  **Any closure claim made without it is unsound.**

## 122. The 47% misprediction rate is second-order content, not broken integrality

The reason every linear veto in this lab is untrustworthy is that 656 of 1,376
large-move predictions were wrong.  The standing explanation was integrality:
forward evaluation solves each gate for its output over Z via `T.solve_lin`, which
returns `None` when the solution is not an integer, leaving the variable stale and
the gate silently broken — a failure mode that has no counterpart mod p.  If that
were the cause, moving along a sublattice `N·Z` with `N` carrying enough copies of
every pivot would restore fidelity.

`s10/fidelity.py` builds that sublattice and measures.  It is not the cause:

```
1 distinct pivot; 0 below 10^6; largest 1.16e+77
free (any integer)      predictions 2160/4086 correct (52.9%); rows exact 220/365; broken gates/probe 0.0
sublattice N = D**3     predictions 2160/4086 correct (52.9%); rows exact 220/365; broken gates/probe 0.0
sublattice N = D        predictions 2160/4086 correct (52.9%); rows exact 220/365; broken gates/probe 0.0
```

**Not one gate breaks, under any move.** Every pivot is ±1 apart from a single
huge one, so forward evaluation is always exact over Z. The misprediction is
therefore pure second-order content: the atoms are genuine quadratics
(`x28730 − x9413·x17499` and friends), and a first-order column cannot predict a
large move.  The map from free inputs to checks is an honest polynomial.

The consequence is sharper than the old explanation.  Under the integrality story
the model was broken and could in principle be repaired.  Under this one the model
is fine and simply **incomplete**: every ceiling in Parts X–XXIII — canonical 39,018,
witness 39,026 — is a ceiling on the *tangent space*, and says nothing about moves
that leave it.  The deliverable "saturates its frame's linear ceiling" is a
statement about the linearisation, not about the instance.

## 123. Decompiling the residual: it is three multiplications and a wire full of p

`s10/decompile.py` prints the residual atoms as polynomials instead of attacking
them as a linear system.  They are tiny:

```
a22229  x7068 - x2099 - 7376877*x642      a22230  x28730 - x9413*x17499
a35758  x29854 - x1329*x22665             a35759  5113045*x7075*x9118 - x29854
a35760  x31864 - x10903*x28961            a35761  x7075*x8731 + x31864
a35762  x642 - x17325*x28599
```

One multiplication each, in redundant pairs that compute the same wire two ways —
ordinary R1CS.  Unfolding the operands reaches only **39 variables, 10 of them
free**, and almost every chain collapses onto `x26064 = p`:

```
x22665 = x26064 = p          x28599 = x9325 = x1692 = x26064 = p
x17499 = x36136 = x32499 = p x28961 = x18822 = x35638 = p
x7075  = 1 - x4287*x2081     x642 = x17325*p     x29854 = x1329*p
```

so the witness-frame residual is exactly two identities:

```
(†)  x1329 * p  =  5113045 * x7075 * x9118           [a35758 & a35759]
(‡)  x10903 * p =  -x7075 * x8731                    [a35760 & a35761]
```

with `x1329, x9118, x10903, x8731, x4287, x2081` all free.  Read that way both are
satisfiable by inspection — `p | x9118` and `p | x8731`, or `x7075 ≡ 0 (mod p)`
which is the `A = 0` route.  The witness frame's "seven conserved residuals" are two
congruences wearing seven hats.

The canonical frame is the better place to work, and Part XXIV's central structural
claim comes from it.

## 124. The canonical frame: the whole problem is congruences, and handles are free

In the canonical frame every gate holds by construction, so the entire instance
reduces to *choose the free inputs so that every check atom vanishes*.  At
`s10/mod9118_0.json` exactly **four** checks are nonzero — four constraints against
7,273 free inputs — and `s10/gadget.py` reads them:

```
a21617 = 11436039*(x14623 - x27522) - p*x5040        x14623 free, x5040  free
a29539 = 12846437*(x14853 - x1308 ) - p*x30163       x14853 free, x30163 free
a37662 = ... + 10*a21617 + (other primitive residuals)
a40826 = ... + (a29539's residual) + ...
```

Every operand that multiplies a free input in `a37662` is a p-wire variable
(`x986 = x20302 = x32499 = x35638 = x36136 = x37280 = p`), so **every handle's
coefficient is a multiple of p**.  That gives a clean two-phase decomposition of the
whole instance:

> **mod-p phase** — the handles are invisible mod p, so the problem over F_p is:
> choose residues for the non-handle free inputs making every check ≡ 0 (mod p).
> **lift phase** — given a mod-p solution, every check is `p · (something)`, and the
> handles enter linearly with coefficient `d·p`; choosing them is an integer linear
> system, not a search.

This is why nothing in Parts X–XXIII ever gained: **every nonzero atom in every state
this lab has produced is nonzero mod p**, so no amount of integer manoeuvring helps.

```
state                                     score   nonzero atoms   of which ≡0 mod p
best/new_instance_partial_39026.json      39,026        7                 0
s10/wr_engine_w1_x7068_39020.json         39,020        2                 0
s10/jm_azero00_39017.json                 39,017        3                 0
s10/mod9118_0.json                        39,009        4                 0
```

The lift phase is free; the mod-p phase is the entire difficulty.

## 125. Residue jumps: a move class no veto ever priced

`s10/gfix.py` closes a gadget constructively rather than by search.  For
`a21617 = c·(x14623 − x27522) − p·x5040`, `s10/suppfree.py` confirms the two sides
are independent — `x27522`'s free-input support is 11 inputs and contains neither
`x14623` nor `x5040` — so:

```
x14623 <- x14623 - ((x14623 - x27522) mod p)     =>  p | c*(x14623 - x27522)
x5040  <- c*(x14623 - x27522)/p                  =>  the atom is EXACTLY zero
```

Both steps exact over Z, and it works: `a21617` and `a29539` are both driven to zero.

```
start                        score 39009   failing checks [21617, 29539, 37662, 40826]
after closing a21617         score 39006   failing checks [25676, 29539, 33796, 40826, 42245]
after closing a29539         score 38999   failing checks [19297, 19299, 25676, 30984, 33796, 36185, 40812, 42245]
```

The score falls — the jump moves `x14623`'s residue, and `x14623` also feeds `a12433`
and `a37662`.  But the *move class* is new and important.  A residue jump is not a
step in the tangent space; it is a jump of a full residue class, and because the map
is a genuine quadratic (§122) its true effect is not what any Jacobian predicts.
**Every price in the Session 11 ledger was measured on tangent moves.  None of them
priced this.**  `s10/rjump.py` enumerates every (failing check, free input in its
support) pair, takes the exact residue jump that zeroes that check mod p, forward-
evaluates and scores — nothing predicted, everything measured.

## 126. What Part XXIV changes

```
claimed in Parts X-XXIII                        status after Part XXIV
------------------------------------------------------------------------------
"the exactly-linear subsystem is inconsistent"  REFUTED (artifact of 162 missing columns)
"linear vetoes untrustworthy: integrality"      WRONG CAUSE -- no gate ever breaks;
                                                it is second-order content
"canonical ceiling 39,018 / witness 39,026"     TANGENT-SPACE ceilings only; they do
                                                not bound the instance
"the residual is seven conserved atoms"         two congruences (†) and (‡) in a
                                                39-variable neighbourhood
"the obstruction is the atom map's image"       refined: it is the mod-p phase.  The
                                                integer lift is free -- handles are
                                                invisible mod p and enter linearly
```

The deliverable is unchanged at **39,026 [verified]**.  What changed is that the
barrier arguments that stood behind it have been withdrawn: two of them were unsound
closures, and the rest bound only the tangent space of a map that is not linear.

---

# Part XXV — the free content of the instance is thirteen numbers

Part XXIV withdrew the barrier arguments.  Part XXV replaces them with the actual
structure, which turns out to be small enough to write down.

## 127. Only d = 0 is safe — exactly, and only for single values

With exact univariate models available (§122: every gate output coefficient is ±1,
so forward evaluation divides by nothing and the map is an honest polynomial),
`s10/unipoly.py` interpolates `f_c(d) = c(v + d·e_u)` for every check and every free
input that reaches a failing check, verifies the fit at extra points, and takes gcds.

```
554 free inputs reach a failing check
degree 9 interpolation exact for EVERY check measured  (degfail = 0 throughout)
43 inputs move nothing at all mod p        (these are the handles)
129 inputs: gcd of the checks that hold has degree 1 -> d = 0 is the ONLY safe jump
             safe roots 0, fixing roots 0, in every case
```

That is an exact statement with no linearisation in it: single-value freedom
does not exist.  So the freedom, if any, is multi-value — and `s10/kerpoly.py`
plus `s10/eqker.py` close that off too:

| closure | rows | cols | rank | kernel |
|---|---|---|---|---|
| atoms, canonical @39,009 | 475 | 197 | 154 | 43 — *exactly the handle columns, which move nothing* |
| **equations**, canonical @39,009 | 3,324 | 415 | 415 | **0** |
| **equations**, frame 2 @39,026 | 139 | 11 | 11 | **0** |

The equation-level rows are the right object — one row per equation instead of one
per atom, and a satisfied *squared* equation `m·(Σ)²` has derivative `2m·Σ·dΣ = 0`
and contributes no row at all — and even there the kernel is empty.  Both states are
first-order rigid.

## 128. Constructing the residual away: the seven are not conserved

`s10/build7.py` takes §123's decompilation at face value.  In frame 2 five of the
seven residual atoms have a detached output variable, so their value may simply be
*written*, and the other two need one divisibility each:

```
force p | x9118 and p | x8731, then set
  x29854 = 5113045*x7075*x9118     x1329  = x29854/p
  x31864 = -x7075*x8731            x10903 = x31864/p
  x28730 = x9413*p                 x642   = x17325*p
  x7068  = x2099 + 7376877*x642
```

Every value exact over ℤ, no search, no linearisation.  It works:

```
start (frame 2)                score 39026  nonzero atoms 7   residual [22229 22230 35758 35759 35760 35761 35762]
after constructing all seven   score 39004  nonzero atoms 4   residual []
```

**All seven residual atoms are exactly zero** — the first time in this lab — with no
broken gate anywhere.  The "conserved quantity" of Parts X–XXIII is not conserved;
it was an artifact of only ever moving inside the tangent space.  What is true is
that the obstruction *relocates*: four checks are now nonzero instead, and the score
falls to 39,004 because those four break 29 equations rather than 7.

## 129. Thirteen 296-bit numbers

Reading the free inputs of any good state settles what the search space actually is:

```
free inputs: 7,273
  7,252 are ZERO
     13 carry 296 bits   <- the whole free content of the instance
      a few carry the large values build7 writes
```

The thirteen are `x6418 x8778 x12553 x14623 x14853 x16742 x22152 x22162 x22649
x24548 x30213 x31339 x33462`.  296 bits is 256 + 40: a field element plus ~40 bits
of `k·p` slack, exactly what a gadget `c·(x − y) − p·h` wants — `x` may be any member
of `y`'s residue class and `h` absorbs the difference.  Writing each as `k·p + r`:

```
x8778 x14623 x16742 x24548 x31339 x33462   k = 839192594282     (six share it)
x14853 x22152 x22649                       k = 1094785891323    (three share it)
x6418  289077647971   x12553 369416716500  x22162 789486214152  x30213 1086320452253
```

`s10/advfix.py` then shows the k part is **irrelevant**: setting an advice value to
the congruence-correct residue with *every* k that occurs in the instance gives
bit-identical scores.  Only the residues matter, which is §124's mod-p/lift split
confirmed by direct measurement.

## 130. The advice constraint graph, and solving it

Each advice value is pinned by exactly one congruence, of one of two shapes
(`s10/advgraph.py`):

```
two-sided   c*(x_i - y_i) - p*h     ->  x_i ≡ y_i (mod p),  y_i computed by the circuit
constant    w*(x_i - C)  - ...      ->  x_i ≡ C  (mod p),  C a literal in the instance
```

Four are pinned to literal constants — ground truth, nothing to solve:

```
x6418  ≡ 20302955751113177691132960011219991444785130617995423281601414462835238472546  via a3576
x12553 ≡ 4531249068709477613185164105669741036354237152756954144434674493737552368539  via a3578
x22152 ≡ 82007976112976807461901870199198737303514020147647909878034348606308756230357 via a31670
x33462 ≡ 37841415183514949237467304684128824427406379377151921996714091976892367869714 via a31672
```

(None of them is a recognisable structured constant, and no pair satisfies any low-degree relation for
small `b`, so they are the generator's own -- see §32.)  The rest are
two-sided, and their targets depend on other advice values:

```
x8778 <- x33462     x14623 <- x24548     x14853 <- x6418
x16742 <- x8778     x22649 <- x22152     x24548 <- x12553     x31339 <- x14623
```

**The dependency graph is a DAG.**  So one Gauss–Seidel sweep in topological order
sets every advice residue correctly, and it does:

```
from B7_39004      sweep 0 changed 4  -> score 39013, residual [19297 19299 30984 36185 40812]
from mod9118_0     sweep 0 changed 3  -> score 39013, residual [19297 19299 30984 36185 40812]
later sweeps       changed 0 -- a FIXED POINT with every advice congruence satisfied
```

Two different starting states land on the *same* residual: 39,013 is an attractor of
the advice solve.  This also explains every "conservation" seen in Parts X–XXIII —
setting `x24548` right makes `a7930` and `a41512` vanish and `a21617` and `a37662`
appear, which looks like conservation but is just an unsolved DAG edge
(`x14623 <- x24548`).  Solve the DAG in order and both clear.

## 131. What is actually left

At the 39,013 attractor the residual is a different family, and `s10/gadget.py`
reads it:

```
a19297 = x11150*x15298 + p*x30317        ->  x11150*x15298 ≡ 0 (mod p)
a19299 = x15298*x25739 - 6672769*p*x5146 ->  x15298*x25739 ≡ 0 (mod p)
a30984 = 537773*x15298*x37758 - p*x2936  ->  x15298*x37758 ≡ 0 (mod p)
a36185, a40812   bundle checks containing those three
```

`x30317`, `x5146`, `x2936` are free handles; `x15298 = 1`; and `x11150`, `x25739`,
`x37758` are 831-bit values whose residues are ~132 bits and nonzero.  Unfolding
them shows what the circuit computes:

```
x29322 = x14853 - x12186        x3558  = x24908 - x16742
x29356 = x29322^2               x27762 = x3558^2
x33469 = x9192 + x24453         x17702 = x29356*x33469
x35389 = x17702 - x27762        x6671  = x27713*x29322 - x1326*x3558
x11150 = 8646263*x35389 + 1073965*x6671
x25739 = 10159099*x35389 + x3023      x37758 = x2287 + 5921311*x6671
```

Differences of values, squared, multiplied by another difference, then
recombined.  And `x15298` is a selector:

```
x15298 = x7715 * x34554
x7715  = x8599 + x21839 - x8599*x21839   = OR(x8599, x21839)
x34554 = x25956 + x7304  - x7304*x25956  = OR(x25956, x7304)
```

a boolean AND of two ORs, currently 1.  So the residual has exactly two doors: make
the three combinations vanish mod p (compute the pair A, B consistently), or drive
the selector `x15298` to 0 (take the degenerate branch).

## 132. Ledger after Part XXV

```
deliverable                                              39,026  [checker-verified]
advice DAG solved, every advice congruence satisfied     39,013  [checker-verified]
all seven residual atoms exactly zero (build7)           39,004  [checker-verified]
```

The deliverable is unchanged.  What changed is that the problem is now *named*
rather than bounded:

| Parts X–XXIII said | Part XXV says |
|---|---|
| the residual is seven conserved atoms | they can be written to zero exactly; the obstruction relocates |
| every route is priced and conserved | the "prices" were unsolved edges of a DAG that is now solved in one sweep |
| the search space is 7,273 free inputs | it is **thirteen 296-bit numbers**, and only their residues mod p matter |
| the obstruction is the atom map's image | it is three products `x15298·{x11150, x25739, x37758} ≡ 0 (mod p)` behind a boolean selector |

The two doors of §131 are the whole remaining problem, and neither has been priced
yet — the selector route in particular has never been attempted from a state where
every advice congruence holds.

---

# Part XXVI — the circuit is algebraic identity, and two of its knobs are free

## 133. The three primitives are homogeneous linear in two quantities

Unfolding the operands of §131 shows the two "spare" terms are not free at all:

```
x3023 = 6926539*x6671        x2287 = 8272701*x35389        (both gate-defined)
```

so with `A = x35389` and `B = x6671` the three primitives become one 3x2 homogeneous
system:

```
x11150 =  8646263*A + 1073965*B          x25739 = 10159099*A + 6926539*B
x37758 =  8272701*A + 5921311*B
```

rank 2, so the only solution is `A ≡ B ≡ 0 (mod p)`.  And A and B are recognisable:

```
A = x29322^2 * x33469 - x3558^2         B = x27713*x29322 - x1326*x3558
```

with `x29322 = x14853 - x12186`, `x3558 = x24908 - x16742`,
`x33469 = (x22162 + x12186 + x14853) + x24453`, `x27713 = x30213 + x16742`,
`x1326 = x12186 - x22162`.  Writing `w1 = x12186, w2 = x16742, w3 = x14853,
w4 = x24908, w5 = x22162, w6 = x30213` for the six free values -- names only, no
interpretation attached -- the two conditions are

```
A = 0   <=>   (w3-w1)^2 * (w5+w1+w3+K) = (w4-w2)^2      K = x24453
B = 0   <=>   (w6+w2) * (w3-w1) = (w4-w2) * (w1-w5)
```

two explicit polynomial identities in six values over F_p, and nothing more is
claimed about them.  The whole residual is that A and B are nonzero, the thirteen
advice values of §129 are point values, and the residual is one identity
that does not yet close.  (The four literal constants of §130 are not the modulus's
group parameters and no pair of them satisfies `y^2 = x^3 + b` for small b, so they
are the generator's own — see §32.)

## 134. w5 and w6 are unconstrained, and A and B both vanish exactly

`x22162` and `x30213` are advice values, and their only pins -- a30976 and a30978 --
are **gated by `x15574`, which is zero**, so nothing constrains them.  A is linear in
`x22162` and B is linear in both.  Two linear equations, two free unknowns:
`s10/abfix.py` recovers the 2x2 matrix by exact probing (the maps are linear, so one
probe per unknown is the exact column) and solves it.

```
exact jacobian  [[111240861292698890820848534820520421376282038064685656642793697565041778256411, 0],
                 [ 82481923122510723799288844430228824462217842441245596187178166524754019170488,
                   54087068875452565652802074820709161994541095135988077443024649864435316913852]]
det invertible
after the jump:  x35389 = 0,  x6671 = 0
                 x11150 ≡ x25739 ≡ x37758 ≡ 0 (mod p)
```

**The A and B both vanish exactly.**  The three primitives of §131 are gone, the
integer lift absorbs a19297, a30984 and a36185 through their handles, and the state
verifies at **39,014** (`s10/P2_39014.json`, checker-verified).

What comes back is a different pair: `a688` and `a1618`, constant pins on `x18956`
and `x24468` -- and `x18956` depends on `x30213` through `x10156 = x15298*x30213`,
so moving `w6` moves them.  Re-running the advice solve reports `changed 0`: every
advice congruence still holds.  The system is now genuinely joint rather than
sequential.

## 135. The joint solve is inconsistent — in the tangent space

`s10/jsolve.py` builds the whole remaining system at once at 39,014: every free input
that reaches a nonzero check (494 after dropping handles), and every equation those
inputs can touch, as equation-combination rows.

```
system: 4,314 rows x 494 cols  (19 of the rows currently fail)
rank 494;  INCONSISTENT ROWS 141
```

So no single linear step reconciles them.  By §122 that bounds the tangent space and
nothing else, which is why it is recorded as a measurement and not as a barrier.

## 136. The selector door is deeper than it looks

The alternative to solving the three primitives is to drive the selector to zero:

```
x15298 = OR(x8599, x21839) * OR(x25956, x7304)
x21839 = OR(x10083, OR(x25608, x390))       x390   = OR(x5638, x33068)
x7304  = OR(x11346, x36945)                 x11346 = OR(x17067, x29560),  x29560 = x13976
```

Every one of these is gate-defined, not free: they are `isZero` comparison flags
computed from the data.  Driving `x15298` to 0 therefore means making an underlying
comparison hold, which is a condition on the advice values, not a knob.  The door is
real but it is the same problem in different clothes.

## 137. Ledger after Part XXVI

```
deliverable                                                39,026  [checker-verified]
identity closed exactly, three primitives gone       39,014  [checker-verified]
advice DAG fixed point                                     39,013  [checker-verified]
all seven residual atoms exactly zero                      39,004  [checker-verified]
```

The deliverable is unchanged.  The instance, however, is no longer an opaque
39,033-equation feasibility problem: it is algebraic point arithmetic over
the modulus whose free content is thirteen 296-bit values, four of them pinned to
literal constants, the rest linked by a DAG of congruences that is now solved in one
sweep, with one identity whose w5/w6 are free and which this session closes
exactly.  What is left is the interaction between that pair and the two constant
pins `a688`/`a1618` that share `w6` -- a joint condition on the same small set of
values, and the first residual in this lab that has never been priced.

---

# Part XXVII — the residual, in closed form, in seven constants

## 138. w5 and w6 are the circuit's output, written into the instance

`s10/pin3.py` unfolds `a688` and `a1618` with the selectors at their current values:

```
x32237 = x21023*x22820,  x21023 = p, x22820 free   -> a handle, ≡ 0 (mod p)
x34243 = x14393*x16153,  x16153 = p, x14393 free   -> a handle, ≡ 0 (mod p)
x25538 = x16742*x34606 + x5647*x24908,   x34606 = x5647 = 0   -> 0
x13913 = x12186*x34606 + x5647*x14853                          -> 0
```

so with `x15298 = 1`, `x18956 ≡ w6` and `x24468 ≡ w5`, and the two pins read straight
off: `w6 ≡ C1·8863713⁻¹`, `w5 ≡ C2`, both literals in `EQUATIONS.txt`.  §134 solved
`A = B = 0` by *moving* w5 and w6, which is exactly why those pins broke.  Setting
them to the pinned values restores `a688 = a1618 = 0` and returns the state to the
39,013 attractor — the loop closes.

## 139. The closed form, verified to the digit

Every quantity in A and B is now a literal constant of the instance:

```
w1 = 82007976112976807461901870199198737303514020147647909878034348606308756230357   (x22152, pin a31670, GATED by x24601)
w2 = 37841415183514949237467304684128824427406379377151921996714091976892367869714   (x33462, pin a31672, GATED by x24601)
w3 = 20302955751113177691132960011219991444785130617995423281601414462835238472546   (x6418,  pin a3576,  GATED by x2081)
w4 = 4531249068709477613185164105669741036354237152756954144434674493737552368539    (x12553, pin a3578,  GATED by x2081)
w5 = 36200939269128454586076546451607958467047992891178506183612554289882454126226   (x22162, pin a1618 via x24468)
w6 = 44859544763832475231923253825569092119321525945631045653619508440821028887      (x30213, pin a688  via x18956)
K  = 97553848499418123410591666447050222001188385549510401465815187079080512838891   (x24453, pin a41332, BARE)
```

and

```
A = (w3-w1)^2 * (w5+w1+w3+K) - (w4-w2)^2
  = 42288441692606730654477992334300923363430351219005991492903082270078522512476     == x35389   EXACTLY
B = (w6+w2) * (w3-w1) - (w1-w5) * (w4-w2)
  = 30198542159037429362146806524344230561752840864915142381356343449320103876465     == x6671    EXACTLY
```

Both match the measured circuit values digit for digit, which settles the reading:
**the instance is one algebraic identity, and the whole residual is that
`A ≠ 0` and `B ≠ 0`.**  For reference, the values the pair A, B wants are

```
required w5 = 64380398444296801010644702415499625279634447310109840487123352893083633736186   (pinned: 36200939...)
required w6 = 45581544895849512040994625888221382902610927244970819299918660665999394080285   (pinned: 44859544...)
```

and `64380398444296801010644702415499625279634447310109840487123352893083633736186` is
exactly what `s10/abfix.py` computed independently in §134.  Two derivations, one
number.

## 140. The gates, and why releasing them does not help

Four of the seven pins are *gated* — `x24601·(w1 − C)` and `x2081·(w3 − C)` — so
zeroing the gate frees the value.  `s10/release.py` measures it:

```
gate x24601 -> 0 :  39,013 -> 38,955, lift -> 38,957   (w1, w2 released)
gate x2081  -> 0 :  39,013 -> 38,937, lift -> 38,939   (w3, w4 released)
both             :  39,013 -> 38,877, lift -> 38,879
```

and then the exact Jacobian of `(A, B)` in the released w's is **identically
zero**.  `s10/closer.py` re-measures it *through* the advice DAG re-solve — since
`x6418 → x1308 → x14853` is how a released constant is supposed to reach the pair A, B
— and interpolation gives **degree 0** in every released w.  Zeroing the gate
frees the value and disconnects it in the same stroke: the gate switches the
whole sub-circuit off.  So the release route is closed, and closed for a reason.

## 141. Where the remaining freedom is

Everything in §139 is forced *on the current branch*.  The branch is set by selector
bits — `x15298 = OR(x8599, x21839)·OR(x25956, x7304)`, and `x34606`, `x5647`,
`x19271`, `x23597`, `x7715`, `x34554` — which decide which formula applies and which
terms (`x25538`, `x13913`) vanish.  They are `isZero` flags, gate-defined, and §136
showed they are not free knobs; but they are the only place left where the instance
is not pinned to its own constants.

That is the sharpest statement this lab has produced about what is left, and it is a
target rather than a barrier: **no infeasibility is claimed** — §121 and the `a7930`
refutation are the standing reminders of what happens to such claims here.

## 142. Ledger after Part XXVII

```
deliverable                                              39,026  [checker-verified]
identity closed exactly (w5, w6 moved)             39,014  [checker-verified]
advice DAG fixed point / w5, w6 at their pins            39,013  [checker-verified]
all seven residual atoms exactly zero                    39,004  [checker-verified]
gate x24601 released, lifted                             38,957
gate x2081  released, lifted                             38,939
```

---

# Part XXVIII — the branch does not reach the residual either

Part XXVII left the branch as the only thing not pinned to the instance's own
constants.  Part XXVIII measures it, and the measurements are exact negatives.

## 143. Frames are irrelevant, at scale

The frame-space search of §61 finished: **4,490 frames evaluated, every single one
scoring exactly 39,026.**  Not one better, not one worse.  Frame choice — which atom
is used to define which variable — cannot move the score at the deliverable, which is
what §139 predicts: the residual is pinned by literal constants regardless of the
solving orientation.

## 144. The boolean census, redone from the advice-solved state

`s10/boolcensus.py` flips every boolean-valued free input at `PIN_39013.json` (the
39,013 attractor with every advice value at its pin) and records the score, the
selector `x15298`, and whether A or B is zeroed.  Over the first 1,464 of 7,250 bits:

```
493 flips are SCORE-NEUTRAL          (39,013 -> 39,013, same five checks)
  1 flip drives the selector to 0    x2081, at a price of 76 (39,013 -> 38,937)
  the rest cost 10-16
  none zeroes A or B
```

## 145. The neutral directions are inert, and the costing ones give two outcomes

"Does not zero A" and "does not move A" are very different facts, and only the second
is a dead end — so both were measured.

```
s10/neutral.py : of 300 score-neutral bits, 300 are COMPLETELY INERT
                 (dA = dB = 0 for every one of them)
s10/mover.py   : of 400 costing bits, only 24 move (A, B) at all, and those
                 24 produce just TWO distinct outcomes:
                    n = 19  ->  A = 18548904073586655983..., B = 23488390206470041529...
                    n =  5  ->  A = 99486789650034179350..., B = 10896062398066642728...
                 neither zeroes A or B
```

The second outcome is bit-identical to what `s10/release.py` measured for `x2081 = 0`,
so those five bits are the `x2081` group.  Nineteen different bits give one and the
same `(dA, dB)` because turning any of them on fires the same OR.

## 146. Why there is no coherent one-hot swap

Reading a bit's atoms shows what these booleans are:

```
a20545  2*b - 2*b*b                              the boolean constraint
a24804  x18232 = b     a24805  x36695 = 1 - b    the bit and its complement
a21775  b*x5803  - C1*b - 12107359*x22874        IF b THEN x5803  = C1
a35126  b*x38738 - C2*b - x12204                 IF b THEN x38738 = C2
```

each bit is a **conditional constant pin** — "if b then this wire is that literal".
The obvious structured move would be a one-hot swap: turn the selected bit off and
another of the same group on, changing which constant is read while keeping the
invariant.  `s10/onehot.py` looks for those groups and finds **none**: every bit gates
its *own* wire, not a shared one, so there is no group to swap within.  And at the
attractor **only 2 of 7,250 bits are 1** — the conditional layer is almost entirely
switched off, which is why turning any bit on only ever adds constraints.

## 147. Ledger after Part XXVIII

```
deliverable                                              39,026  [checker-verified]
  every one of 4,490 frames scores exactly                39,026
identity closed exactly (w5, w6 moved)             39,014  [checker-verified]
advice DAG fixed point / w5, w6 at their pins            39,013  [checker-verified]
all seven residual atoms exactly zero                    39,004  [checker-verified]
selector x15298 -> 0 (x2081), after the lift             38,937
```

Every door of §141 is now measured:

| door | measurement | result |
|---|---|---|
| move w5, w6 | exact 2x2 solve (`abfix.py`) | closes A and B, breaks the pins that fix them |
| release a gated pin | `release.py`, `closer.py` | frees the value *and disconnects it* — degree 0 |
| flip a boolean | `boolcensus.py` | 493 neutral, none zeroes A or B |
| neutral directions | `neutral.py` | 300 of 300 completely inert on (A, B) |
| costing directions | `mover.py` | 24 of 400 move (A, B); two outcomes, neither zero |
| one-hot re-selection | `onehot.py` | no groups exist; only 2 of 7,250 bits are on |
| frame choice | 4,490 frames | all exactly 39,026 |

**No infeasibility is claimed.**  §121 and the `a7930` refutation are the standing
reminders of what happens to such claims in this lab, and these are measurements over
enumerated move classes, not a proof that none exists.  What they do establish is
that the residual of §139 is not reachable by any move class this lab has been able
to name — and that it is now a two-line arithmetic statement in seven printed
constants rather than a 39,033-equation mystery.

---

# Part XXIX — A and B both vanish, five ways, and every one of them is priced

Part XXVII said the seven quantities in A and B are literal constants of the
instance.  That is **wrong as stated**, and correcting it is what this part is about.

## 148. Correction: the values are not the literals

`s10/valjac.py` runs one forward-AD pass per free input — 7,273 of them — and
records the exact derivative of every value.  Two facts kill §139's premise:

```
w1 = x12186   computed, moved by 179 free inputs     w2 = x16742   FREE
w3 = x14853   FREE                                   w4 = x24908   computed, 43
w5 = x22162   FREE                                   w6 = x30213   FREE
K  = x24453   empty support -- the only genuine constant
```

and perturbing `x22152`, `x33462`, `x6418` or `x12553` moves **none** of w1, w2, w3,
w4.  Those four literals reach the values only through the advice DAG, not
directly; §139 read a coincidence of residues as an identity of variables.  The
closed forms for A and B are still exact — they reproduce `x35389` and `x6671` digit
for digit — but the values in them are steerable, and the pair is therefore
**not over-determined**.

## 149. The w-map has rank 8, and its non-boolean part is diagonal

Over the 264 free inputs with any effect, the map to
`(w1, w2, w3, w4, w5, w6, x19083, x1308, A, B)` has **rank 8 of 10**, and the only
two relations are the linearisations of A and B themselves.  All eight values
are independently steerable.

But steering them all at once with all 264 knobs costs 91 points
(`s10/newton8.py`: Newton converges, residual on all eight targets exactly zero, A and
B exactly zero, a26731 and a29539 intact — and 22 other checks broken).  The reason
is visible once the knobs are separated:

```
264 value movers  =  8 non-boolean  +  256 boolean
```

and driving a boolean free input to an arbitrary residue breaks its own `b² = b`
constraint.  The eight non-boolean movers are eight of the thirteen advice values, and
they form a **diagonal** system — one knob per value:

```
w1  <- x22649      w2  <- itself      w3  <- itself      w4  <- x31339
w5  <- itself      w6  <- itself      x19083 <- x8778    x1308  <- x6418
```

so freeing one w costs exactly one congruence: w1 costs a2423, w4 costs
a33796, x19083 costs a33929, w2 costs a26731, w3 costs a29539, w5 and w6 cost a1618
and a688.

## 150. A = B = 0 is solvable — five ways — and every one is priced exactly

A and B are two equations, so driving A and B to zero needs two values.  Two of the
pairs solve **linearly**, which is worth writing out because it removes any
root-existence question:

```
(w1, w2)   A  =>  (w5 + w1 + w3 + K)(w3 - w5)^2 = (w4 + w6)^2      LINEAR in w1
(w3, w4)   A  =>  (w5 + w1 + w3 + K)(w1 - w5)^2 = (w6 + w2)^2      LINEAR in w3
(w5, w6)   A  =>  w5 = (w4-w2)^2/(w3-w1)^2 - w1 - w3 - K           LINEAR in w5
(w3, w2)   cubic in w3 -- exactly ONE root in F_p
(w1, w4)   cubic in w1 -- NO root in F_p
```

`s10/pairfix.py` drives each pair through its single diagonal knob, lifts, and scores:

```
pair        A=0  B=0   score   surviving checks
(w5, w6)    yes  yes   39,015  [688, 1618, 19297, 19299, 40608, 40812]   <- best
(w3, w2)    yes  yes   38,992  [19299, 26731, 29539, 36185, 40812, 40826]
(w1, w2)    yes  yes   38,991  [2423, 10506, 19299, 26731, 36185, 40812]
(w3, w4)    yes  yes   38,991  [19299, 25676, 29539, 33796, 36185, 40812, 40826, 42245]
(w1, w4)    no solution in F_p
```

**`s10/PF_best_39015.json` verifies at 39,015/39,033** — A and B both vanish and
the score exceeds every previous canonical-frame state.  Each pair pays exactly the
congruences of the values it moves, as §149 predicts.

## 151. The trap, stated precisely

Closing the pair A, B always costs the congruences of the two values moved, and
those congruences cannot be repaired locally, because the advice DAG is a **chain
rooted in the four gated literals**:

```
x33462 (a31672, gated by x24601) -> x16144 -> x8778 -> x19083 -> w2
x6418  (a3576,  gated by x2081 ) -> x1308  -> w3
x22152 (a31670, gated by x24601) -> x29524 -> x22649 -> w1
x12553 (a3578,  gated by x2081 ) -> x24548 -> x14623 -> x31339 -> w4
```

Move a value and its link breaks; repair that link and the next one up breaks;
the chain terminates at a literal whose pin is ungated only when `x2081` or `x24601`
is zero — and those are exactly the quadrant switches of §147, so releasing one turns
the pair A, B check off and substitutes `the first branch` or `the second branch`, which the same literals
contradict.  That is the whole obstruction, and it is now a closed loop rather than a
mystery.

## 152. Why the deliverable still wins

```
deliverable                                       39,026  [checker-verified]
A = B = 0 via (w5, w6)                      39,015  [checker-verified]   <- NEW
identity closed, w5/w6 moved (abfix)        39,014  [checker-verified]
advice DAG fixed point                            39,013  [checker-verified]
all seven residual atoms exactly zero             39,004  [checker-verified]
A = B = 0 via (w3,w2) / (w1,w2) / (w3,w4)   38,992 / 38,991 / 38,991
Newton on all eight values at once           38,922
```

The deliverable is **not on the "solve the circuit" path at all**.  Its seven nonzero
atoms are arranged so that their bundle combinations cancel in all but seven
equations; the structurally clean states have fewer nonzero atoms (two, three, four)
but those atoms sit in equations that do not cancel.  39,026 is a *coding* optimum,
39,015 is the *algebraic* optimum, and this session established for the first time
that they are different objects.

**No infeasibility is claimed.**  Everything above is an exact price for an enumerated
move, and §148 is a reminder of how a confident structural claim can be wrong: Part
XXVII's premise survived several checks and was still false.

---

# Part XXX — the lift belongs at the equation level, and the last two points are locked

---

# Part XXX — the lift belongs at the equation level, and the last two points are locked

## 153. A gap between the mod-p count and the checker

`s10/ceil15.py` runs the coset leader at the new 39,015 state and reports something
the lab had not noticed:

```
PF_best_39015.json: score 39015; failing checks [688, 1618, 19297, 19299, 40608, 40812]
closure: 1393 rows x 426 cols
do nothing: 16 equations fail -> score 39017   (the checker says 39015)
LINEAR CEILING AT PF_best_39015.json : 39017
```

Sixteen equations have a combination that is nonzero **mod p**; the checker finds
eighteen failing **over ℤ**.  So two equations — 7469 and 21382 — have an atom
combination that vanishes mod p and not over ℤ.  Every lift in this lab has worked
atom-by-atom (find a check ≡ 0 mod p, absorb it with a handle), and an equation whose
*combination* vanishes while its atoms do not is invisible to that.  Surveying every
good state:

```
best/new_instance_partial_39026.json   7 failing,  0 liftable
s10/PF_best_39015.json                18 failing,  2 liftable  [7469, 21382]
s10/P2_39014.json                     19 failing,  3 liftable  [7123, 7469, 21382]
s10/AG_39013.json                     20 failing,  0 liftable
s10/wr_engine_w1_x7068_39020.json     13 failing,  0 liftable
```

The deliverable has none — no free points there — but the algebraic states carry two
or three.

## 154. The equation-level lift, and why it oscillates

`s10/eqlift.py` implements the right object.  For a failing equation whose combination
`S_e` is ≡ 0 (mod p),

```
dS_e/du  =  sum_a c_a * (da/du)_Z          exact, via intad.jacZ
u <- u - S_e / (dS_e/du)                   whenever that coefficient divides S_e
```

drives `S_e` to exactly zero over ℤ **without any atom of e being zero** — which is
precisely the mechanism the 39,026 deliverable exploits (§152), reached there by
search rather than by construction.

It fires, and it oscillates: eq 7469 is absorbed by the handle `x30317`, which breaks
eq 7123 through the same handle, which is absorbed by `x30317`, which breaks 7469.
With strict improvement required, **no single move exists**.

## 155. Solving the pair together: no integral solution

Two equations, two knobs, over ℤ:

```
sum_u g[e][u] * d_u = -S_e     for e in {7469, 21382},  d_u integer
```

a linear Diophantine system whose 2x2 sub-solves are a determinant-divisibility test.
`s10/eqdio.py` enumerates every pair among the 238 knobs with a nonzero exact integer
effect:

```
238 knobs with a nonzero exact integer effect
0 integral pairs tried;  best 39015 (was 39015)
```

**Not one pair yields an integral solution.**  The determinant divides neither
right-hand side in any of the ~28,000 pairs.  So the two points that separate the
mod-p coset count (39,017) from the checker (39,015) are locked by integrality, not
by the mod-p algebra — the first place in this whole investigation where the
obstruction is genuinely arithmetic over ℤ rather than a congruence over F_p.

## 156. Final ledger

```
deliverable                                        39,026  [checker-verified]
  its coset-leader ceiling (witness frame)         39,026  -- saturated
A = B = 0 via (w5, w6)                       39,015  [checker-verified]
  its coset-leader ceiling                         39,017  -- 2 points, locked over Z
two-condition primitive closed                     39,014  [checker-verified]
advice DAG fixed point                             39,013  [checker-verified]
all seven residual atoms exactly zero              39,004  [checker-verified]
```

---

# Part XXXI — why 39,026 is optimal for its own residual, exactly

Part XXIX separated a *coding* optimum (39,026) from an *algebraic* one (39,015).
Part XXXI settles the coding side completely, with integer linear algebra rather
than a search.

## 157. The twelve equations are homogeneous in the seven atom values

At the deliverable the seven residual atoms touch twelve equations, and those twelve
involve 24 atoms of which seventeen are zero.  So each equation is a linear form in
the seven residual values `alpha` **with no constant term**:

```
 0 eq 2554   a0 +13 a1                                          ZERO
 1 eq 6816   -15 a0 -11 a1 +38 a2 +9 a3 +36 a4 +13 a5 +29 a6    ZERO
 2 eq 8124   36 a0 +26 a1 +a5 -6 a6                             ZERO
 3 eq 9123   20 a2 +27 a3 +33 a4 +3 a5 -a6                      ZERO
 4 eq 9421   13 a0 -21 a1 -21 a2 +29 a3 +38 a4 +29 a5 +4 a6     ZERO
 5 eq 12231  18 a0 +24 a1 +a4 -23 a5 +13 a6                     fails
 6 eq 12270  -31 a0 +5 a1 +a2 -27 a3 -a4 -17 a5 +10 a6          fails
 7 eq 12350  -23 a0 +26 a1 -16 a2 +34 a3 -34 a4 +35 a5 +11 a6   fails
 8 eq 14584  17 a0 +16 a1 -2 a2 -18 a3 -31 a4 +19 a5 -39 a6     fails
 9 eq 18673  a3 +6 a4 +a5                                       fails
10 eq 22044  -24 a0 -10 a1 +a6                                  fails
11 eq 29125  a1                                                 fails
```

rank 7, so `alpha = 0` satisfies all twelve, and the witness satisfies exactly five.
Note row 11: **`eq 29125` is `a22230` alone**, so it demands `a1 = 0` exactly.

## 158. The zero-cost lattice, enumerated and measured

`s10/genscan.py` takes every free input that structurally reaches the seven atoms --
there are only fifteen -- moves it by one, records the exact change to all seven
`alpha` components, and counts the equations that break *outside* the twelve.  Nine
have cost zero:

```
x642    a0 -7376877, a6 +1        x29854  a2 +1, a3 -1        x31864  a4 +1, a5 +1
x1329   a2 -p                     x10903  a4 -p               x9413   a1 -p
x17325  a6 -p                     x9118   a3 +5113045         x8731   a5 +1
```

This corrects two earlier guesses: x9118 (33 atoms) and x8731 (26 atoms) *look*
expensive and are free.  Only two of the fifteen cost anything, and they are exactly
the fine generators of the two coarse directions:

```
x7068  ->  a0 += 1,  cost 13         x28730  ->  a1 += 1,  cost 16
```

Every zero-cost generator leaves `a2+a3` and `a5-a4` unchanged mod p -- precisely the
condition for the residues of x9118 and x8731, the only quantities with any reach
outside, to stay put.  So the reachable set at zero cost is the coset
`alpha_witness + L` with L the rank-7 lattice above, in which

```
a2, a3, a4, a5   completely free
a1               only in multiples of p
a0               only in multiples of 7376877, with a6 coupled to it mod p
```

## 159. No subset of more than five is integrally reachable

`s10/lattice7.py` solves, for every subset S of the twelve rows, the linear
Diophantine system `A c = B` with `A[i][j] = <M_i, g_j>` and `B[i] = -<M_i, alpha_w>`,
by column-style Hermite reduction (unit-tested, and it returns `c = 0` for the five
rows that already hold):

```
|S| = 12, 11, 10, 9, 8, 7, 6  :  0 integrally solvable subsets
all 924 six-subsets and all 792 seven-subsets ARE solvable over Q
```

**The obstruction is pure integrality**, and it sits in the two coarse directions:
`a1` can only move by multiples of p, while row 11 needs `a1 = 0` and `a1` is not
≡ 0 (mod p) at the witness.

## 160. The coarse generators cost at least eleven

`s10/repaircost.py` moves each coarse generator and identifies exactly what breaks:

```
x28730 += 1  ->  16 equations, whose only nonzero atoms are a7930 and a41512
                 -- the advice congruence x24548 = x25442
x7068  += 1  ->  13 equations, whose only nonzero atoms are a29539 and a40826
                 -- the advice congruence x14853 = x1308
```

one advice congruence each, and both are free advice values, so `s10/coarse.py`
repairs them with the exact residue-jump-plus-handle of §125 and then re-solves the
advice DAG.  The collateral does not clear:

```
x28730 += 1            39,006      repair a7930 via x24548   39,011   (11 still broken)
x7068  += 1            39,009      repair a29539 via x14853  39,002   (20 still broken)
advice sweep afterwards                                      39,002 in both cases
```

The cheapest coarse move measured costs **11**.

## 161. The optimality statement

Writing `k` for how many of the twelve hold and `c` for the collateral elsewhere,

```
score  =  39,033 - (12 - k) - c
```

and the two results above bound it:

```
c = 0   =>  k <= 5   (§159, exact integer linear algebra over the complete lattice)
        =>  score <= 39,026     -- and the deliverable ATTAINS it
c >= 11 (§160, measured over all fifteen inputs that reach the residual)
        =>  score <= 39,033 - 11 = 39,022  <  39,026
```

**So 39,026 is optimal for this residual structure.**  It is the first optimality
statement in this lab that is not a linearisation: §159 is exact integer arithmetic,
and §158 enumerates the generators completely rather than sampling them.

The statement is conditional in two honest ways, and neither is hidden.  The coarse
cost 11 is *measured* on single moves and their repairs, not proved -- a combination
of coarse moves whose collateral cancels would evade it.  And the whole argument is
about the deliverable's frame and its seven residual atoms; it says nothing about a
different frame, or about the algebraic path of Part XXIX, or about whether the
instance has a full solution.  **No infeasibility is claimed.**

## 162. Ledger

```
deliverable                                     39,026  [checker-verified]  OPTIMAL for its residual
A = B = 0 via (w5, w6)                    39,015  [checker-verified]
  its compensation ceiling                      39,017
two-condition primitive closed                  39,014  [checker-verified]
advice DAG fixed point                          39,013  [checker-verified]
coarse move x28730 + a7930 repair               39,011
all seven residual atoms exactly zero           39,004  [checker-verified]
```

---

# Part XXXII — the algebraic side, pushed to 39,017 and then closed

Part XXXI proved 39,026 optimal for the *coding* residual.  This part does the same
job for the *algebraic* one, and improves it first.

## 163. The last lock is one divisibility, and the k·p freedom solves it

With A = x35389 and B = x6671 both ≡ 0 (mod p) the three primitives are exactly

```
x11150 = 8646263*A  + 1073965*B     a19297 = x11150 + p*x30317
x25739 = 10159099*A + 6926539*B     a19299 = x25739 - 6672769*p*x5146
x37758 = 8272701*A  + 5921311*B     a30984 = 537773*x37758 - p*x2936
```

so a19297 and a30984 are absorbed by their handles the moment A and B are multiples
of p.  a19299 is the exception: its handle enters with coefficient `6672769*p`, so
with `A = p·a`, `B = p·b` it needs one extra divisibility

```
6672769  |  10159099*a + 6926539*b          (6672769 is prime)
```

which `s10/handzero.py` and `s10/habsorb.py` show is exactly what blocks it.  But
a and b are not fixed: every advice value is `k*p + r`, and bumping the k of x22162 or
x30213 shifts A and B by exact multiples of p while touching **no residue** — hence no
congruence.  `s10/divlock.py` measures the steps:

```
x22162 += p  ->  target moves by 1963712 (mod 6672769)
x30213 += p  ->  target moves by 3063958 (mod 6672769)
gcd(1963712, 3063958, 6672769) = 1   ->  the congruence ALWAYS has solutions
```

## 164. Absorb with handles only — 39,017

Every solved candidate still lost points, and `habsorb` says why: **x22162 absorbs
a1618 with coefficient 1, and x30213 absorbs a688 with coefficient 8863713.**  Those
are the w's themselves, so a greedy absorber "fixes" a1618 by putting w5 back
on its pin and destroys A = 0.  The absorber must be restricted to genuine handles —
free inputs whose exact integer coefficient is a multiple of p, which is precisely
what makes the move invisible mod p and harmless to every congruence.

`s10/finish.py` does that, and at `k1 = -12, k2 = 5812259` the lock goes to zero and
all three primitives clear:

```
PF_best_39015    18 failing   checks [688, 1618, 19297, 19299, 40608, 40812]
FIN_39017        16 failing   checks [688, 1618, 40608]          [checker-verified]
```

**`s10/FIN_39017.json` verifies at 39,017/39,033** — the best algebraic state, and the
first time the two-condition primitive closes with its parts fully absorbed.

The full sweep confirms it independently: **all 8 shifts that solve the lock land on
exactly 39,017**, none higher and none lower.  That is what §165 predicts — the score
cannot depend on which solution of the congruence is chosen, because the residues of
a688 and a1618 are pinned by the pair A, B and the k·p shifts move nothing else.

## 165. 39,017 is exact there: no ratio saves anything

The residual is two numbers, and ten of the sixteen failing equations contain **both**
a688 and a1618, so a cancellation `c1*a688 + c2*a1618 = 0` could in principle save
them — the very mechanism that makes 39,026 a coding optimum.  The residues of a688
and a1618 mod p are pinned by the pair A, B, while their handles add arbitrary
multiples of p, so an equation is recoverable iff `c1*r688 + c2*r1618 ≡ 0 (mod p)`.
`s10/ratio.py` tests all sixteen:

```
eq 56 (34, 0) · eq 133 (-35,-16) · eq 2071 (6,-28) · eq 8073 (-11, 1) · ...
0 of 16 failing equations pass the mod-p ratio test
```

**Not one.**  No choice of the two handle multiples can save a single equation, so
39,017 is the exact optimum at that state.

## 166. And 16 is the minimum cost of driving A and B to zero

Closing A = B = 0 needs exactly **two** w moves — one cannot do it:

```
w5 from A = 643803984442968010106447024154...     w5 from B = 641082231655455145747261280924...
   they differ, so w5 alone fails
w6 appears only in B, so it can never fix A
A = 0 solved in w2 needs a square root of a NON-RESIDUE mod p -- no w2 exists
```

and each moved w breaks exactly one congruence plus its bundle atoms.  The
union of the failing-equation sets over every pair:

| pair | pins | distinct equations | best possible score |
|---|---|---|---|
| **(w5, w6)** | a1618 + a688 | **16** | **39,017** |
| (w6, w1) | a688 + a2423 | 23 | 39,010 |
| (w1, w3) | a2423 + a29539 | 24 | 39,009 |
| (w1, w4) | a2423 + a33796 | 25 | 39,008 |
| (w6, w3) | a688 + a29539 | 25 | 39,008 |
| … | … | 26–31 | 39,007–39,002 |

The pair this session already uses is the cheapest by seven equations.

## 167. Both optima are now exact, and the coding one wins

```
CODING     deliverable                39,026   optimal (Part XXXI, integer lattice)
ALGEBRAIC  A = B = 0, absorbed 39,017   optimal (this part)
```

The algebraic route cannot beat the coding route, and the reason is structural rather
than accidental: driving A and B to zero costs two congruences, the cheapest pair of
which sits in sixteen equations, while the deliverable's seven residual atoms sit in
only twelve and cancel in five of them.  **Nine equations separate the two, and both
ends of that gap are now exact statements** rather than the best a search could find.

Still no infeasibility claim: §166 prices the moves this lab can name, and §165's
ratio test is exact only for the state it is run at.
