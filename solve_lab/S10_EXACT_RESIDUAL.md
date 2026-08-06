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
