# The irreducible residual problem — size, shape, and solver suitability

All numbers below are **measured** on the instance, not estimated.
Best verified assignment: `best/new_instance_partial_39024.json` = 39,024/39,033 (9 failing).

## 1. The reduction chain

| stage | size |
|---|---|
| raw instance | 39,033 equations, 38,748 unknowns |
| shared sub-expressions ("atoms") | 42,267 (deg 1: 19,780 · deg 2: 21,788 · deg 4: 699) |
| after recovering the setter's gate orientation | 31,475 gates, **7,273 free inputs**, 10,792 checks |
| free inputs by type | **1,156 boolean**, ~6,117 integer-valued (256–300 bit) |
| residual after forward construction | **2 scalar congruences mod p** |

with `p = 2^256 − 2^32 − 977` (the secp256k1 field prime, which appears *exactly* as the value of
the 220-variable identity wire).

**The whole instance is satisfied iff**

```
D1 := x_7068 − x_2099 − 7376877·x_642   ≡ 0   (mod p)
D2 := x_4432 − x_19964 − x_28730        ≡ 0   (mod p)
```

everything else is satisfied by construction. So the irreducible problem is **512 bits of
constraint** (two 256-bit congruences).

## 2. How much of that is COMBINATORIAL (the annealer question)

This is the decisive measurement. Of the 1,156 boolean free inputs:

* **2** move `(D1, D2)` mod p at all — `x_2081` and `x_4287`, the two MUX control bits;
* their deltas are **identical**, so the rank of the bit → residue map over GF(p) is **1 of 2**;
* the other 1,154 leave both residues **exactly** unchanged (902 of them are wholly inert — they
  change 3 variables and no atom's zero/nonzero status).

Of the ~6,117 integer free inputs, exactly **4** move the residues — `x_6418`, `x_12553`
(the MUX data inputs) and `x_9413`, `x_17325` (quotient handles). In any *complete* solution
`x_6418`/`x_12553` are pinned mod p by their load gates and `x_9413`/`x_17325` contribute only
multiples of p, i.e. nothing mod p.

> **Consequence: the 256-bit message is provably decoupled from the verification residues.
> There is no subset-sum / knapsack / MAX-CUT kernel hiding in this instance.** The reachable set
> of `(D1, D2)` mod p is a handful of points, none of them `(0,0)`.

The hardness is **arithmetic** (hit a prescribed 256-bit residue through a modular circuit), not
combinatorial.

## 3. The smallest self-contained algebraic statement

After eliminating every wire, the kernel is six quantities mod p and two quadratic forms:

```
S = A·u² − w²                        T = B·u − w·c
L1 = 8646263·S + 1073965·T           L2 = 10159099·S + 6926539·T
L3 = 8272701·S + 5921311·T
require   p | L1,  p | L3,  6672769·p | L2      ⟺   S ≡ T ≡ 0 (mod p)  + one condition mod 6672769
```

with `u = x_29322`, `w = x_3558`, `A = x_33469`, `B = x_27713`, `c = x_1326`.
Eliminating `w` gives `u²·(A·c² − B²) ≡ 0`, so **either `u ≡ w ≡ 0`, or `A·c² ≡ B² (mod p)` — and
the second branch is a pure Legendre-symbol condition: it is open iff `A` is a quadratic residue.**

That is the natural place for number theory, and it is cheap to test (one Legendre symbol per
configuration). Measured: at the *forced* configuration (`x_14853` pinned to `K1`) `A` is a
**non-residue** for both reachable values of `x_12186`, which closes that branch. A ~50% sample of
hypothetical shifts of `x_14853` *would* be residues — the setter's pin is precisely what closes it.

## 4. Why quantum annealers are the wrong machine here

Two independent blockers, the second fatal:

**(a) Size.** A QUBO for this must encode 256-bit modular arithmetic. One 256-bit modular multiply
costs ≈ 256² ≈ 65,000 binary variables for the partial-product array plus carry chains; the residual
kernel above contains at least 4 such multiplies (`A·u²`, `w²`, `B·u`, `w·c`) and the chains feeding
`A, B, c, u, w` add more. Realistic logical size: **10⁵–10⁶ binary variables** before
minor-embedding. D-Wave Advantage2 offers ~4,400 qubits at ~20-way connectivity, and embedding a
dense arithmetic circuit costs a further quadratic blow-up — so ~10⁸–10¹⁰ physical qubits.

**(b) Precision — the real blocker.** Even granting unlimited qubits, the objective
`(Σ aᵢbᵢ − t − k·p)²` with 256-bit coefficients needs roughly **512 bits of dynamic range in the
coupler weights**. Current annealers deliver ~4–6 *effective* bits after analog control error
(~1–2% ICE). No amount of qubit growth fixes this; the energy landscape simply cannot be
represented. This is the same reason annealers have not threatened 256-bit factoring or discrete
log despite repeated attempts (published D-Wave factoring records are ~6-digit numbers, and those
use heavy problem-specific preprocessing).

**(c) There is no low-precision residue to peel off.** If the message bits *had* controlled the
residues, one could hand a 256-variable subset-sum to an annealer. Measurement (§2) rules that out:
rank 1 of 2, from the two MUX bits.

## 5. What *is* small enough to hand to a solver

| object | size | status |
|---|---|---|
| MUX quadrant choice | 4 configurations | enumerated exhaustively |
| knob-lattice integer program (what produced 39,024) | 13 equations × 8 integer knobs, 4 mod-p congruences | solved exactly (SNF); optimum 9 failing |
| single-bit / two-bit message search | 1,156 and 667,590 states | exhausted, no improvement |
| the residual kernel | 6 unknowns mod p, 2 quadratic forms | **open** |
| the closed-form identity of `S9_STRUCTURE.md` §12 | 1 scalar equality mod p | **open — this is the target** |

## 6. Recommended direction

Not annealing. The tractable shapes are:

1. **Legendre / quadratic-residue conditions** — the second core branch opens iff `A` is a QR; find
   a reachable configuration where the symbol flips. Cost: one modular exponentiation per test.
2. **Lattice reduction (CVP/LLL)** on the accessible knobs *if* the map to the residue is affine
   mod p over a large enough knob set — this was never run, and is the standard tool at this size.
3. **The §12 identity directly**: make the residue that `x_24601`'s pin puts on `x_22649` agree with
   the residue the MUX-branch obligation forces on `x_9118`. Both sides are explicit; the question
   is whether any reachable configuration makes them coincide.
