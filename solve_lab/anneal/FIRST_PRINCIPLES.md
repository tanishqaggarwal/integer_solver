# A number-theoretic reading, and how small the encoding can actually get

This is the second pass. The first (`ENCODING.md`) built a faithful QUBO and
measured it. This one asks the sharper question: *what is the encoding actually
paying for, and is there a reading of the arithmetic that avoids the bill?*

---

## 1. The instance has no arithmetic weakness (`weakness.py`)

Every standard way a curve is secretly easy, checked rather than assumed:

| test | result |
|---|---|
| group order `n` prime | **yes** — no Pohlig–Hellman, no small subgroups |
| anomalous (`n = p`, Smart's attack, polynomial time) | no |
| MOV/Frey–Rück embedding degree `e` with `p^e ≡ 1 (mod n)` | **> 200** — the pairing lands in a >51,200-bit field |
| `j = 0`, CM by `Z[ζ₃]`, order-6 automorphism group | present — buys `√6 ≈ 2.4×` in Pollard rho, i.e. `2^128 → 2^126.7` |
| base point `G` = standard secp256k1 generator | no — `G` is randomised too |
| `T = c·G` for small `c` | nothing under 64 |

So nothing here is broken. The remaining hope was structure in *how the instance
was generated*, which is the subject of the parallel search tracks.

## 2. What the encoding is really paying for

`E(F_p)` is cyclic of prime order `n`, so **as a group it simply is `Z_n`**. Under
that isomorphism the entire instance is a modular subset-sum:

```
find c_j ∈ {0,1}  with  Σ c_j·ℓ_j ≡ ℓ_T  (mod n),     ℓ_j = dlog_G(Q_j)
```

which needs no field arithmetic at all. The isomorphism `E → Z_n` is exactly the
discrete logarithm — so the split is clean and total:

> **Everything above a modular subset-sum is the price of not knowing the logs.**

This is worth stating because it says precisely where the cost lives. We *do* know
256 logs already — `ℓ(P_i) = 2^i` is what the doubling chain means — but the one
log we need is `ℓ_T = k`, the answer. Any additional relations with known
logarithms would collapse the encoding. For prime-field curves, producing such
relations is the open problem (there is no index calculus); the encoding is not.

## 3. Both halves of that split, measured (`subsetsum.py`, `embed.py`)

The binding hardware constraint is not the logical variable count — it is the
largest **clique** in the Hamiltonian. Every penalty here is a square of a linear
form, and a square over `c` terms makes a `K_c`. On Pegasus, `K_182` consumes all
5,760 qubits, so a `K_c` costs ≈ `c/6` physical qubits per logical variable.

**If the logs were known** — modular subset-sum over 256 relations (256 binary
variables are information-theoretically required: `2^m` subsets must cover `2^256`
targets, so `m < 256` has no solution at all):

| variant | logical | max clique | ≈ physical |
|---|---|---|---|
| no compression | 2,371 | 169 | 66,783 |
| chunk-compressed | 6,376 | 71 | 75,449 |
| full-adder compressed | 67,739 | 5 | 67,739 |

The product is **invariant at ≈ 6.7 × 10⁴ physical qubits** — shrinking the clique
inflates the ancillas by the same factor. That is ~15× a current annealer: within
one hardware generation.

**As the instance actually stands** — the full comb ladder:

| encoding | logical | max clique | ≈ physical | `\|J\|` |
|---|---|---|---|---|
| binary carries, `w = 9` | 9.08 × 10⁶ | 139 | 2.1 × 10⁸ | 2²¹ |
| **full-adder carries, `w = 8`** | 2.66 × 10⁷ | **5** | **2.66 × 10⁷** | **2⁹** |

**This reverses the ranking in `ENCODING.md`.** Binary carries use 3× fewer
logical qubits but 8× more physical ones, because their column equations make
`K_139`. Measured in the currency that matters, the full-adder encoding is the
minimum: **2.66 × 10⁷ physical qubits, ~6,000× current hardware.**

One concrete win got us there: the "exactly one digit" constraint was
`(Σu − 1)²`, a `K_D` with `D = 2^w` — that single penalty was dominating the whole
embedding. Rewritten as a sequential-counter prefix chain (`p_t = p_{t-1} + u_t`,
3-term penalties) the largest clique in the entire Hamiltonian drops to **5**.
Faithfulness re-verified exhaustively after the change.

## 4. Why the comb is already optimal (`optimality.py`)

The measured per-window cost obeys a clean two-parameter law, accurate to **0.8%**
across `w = 1…11`:

```
per_window(w) = 268,606 + 83.0·2^w        (binary)
              = 698,240 + 513.2·2^w       (full-adder)
```

giving one 256-bit modular multiplication ≈ **89,535** qubits and one table entry
≈ 83.

*Coverage bound.* A comb with `A` additions and tables of size `D` offers `D^(A+1)`
digit tuples, so addressing every `k < 2^256` forces `A ≥ 256/log₂D − 1 = 256/w − 1`.
Additions cannot be traded away faster than logarithmically while look-up cost
grows as `2^w` — which is exactly the measured trade-off, and why the optimum sits
at `w = 8–9` rather than at either extreme. Separately, `A ≈ 29` additions is at
the Pippenger bound for a 256-bit scalar with unbounded precomputation, so the
group-operation count is not improvable either.

*The floor.* Three multiplications per addition is minimal for the affine law
(`λ·d = e`, `λ² = x₃+x₁+x₂`, `λ·(x₁−x₃) = y₃+y₁`; projective coordinates trade the
inverse for *more* products). So within this family

```
cost(w) ≥ (256/w)·3M,   M = 89,535 measured
```

and getting that under 4,400 needs `w ≥ 15,628` — a table of `2^15628` points. The
floor is **one modular multiplication, ≈ 9 × 10⁴ qubits, ≈ 20× a real annealer** —
and no encoding of this decision problem contains fewer than one, because deciding
it means comparing `x(kG)` with `x(T)` in `F_p`.

## 5. Where this leaves the multi-run question

The user's constraint was: many runs are fine provided each fits the hardware.
Section 4 says the smallest *sound* indivisible unit is one modular
multiplication at ~9 × 10⁴ qubits. Below that you are no longer decomposing the
problem, you are clamping variables and hoping — and `ENCODING.md` §6 measured
that hope failing: simulated annealing cannot settle a single modular
multiplication past 8-bit words, nor fill uniquely-determined ancillas on a 4-bit
ladder at 100,000 sweeps.

So the honest summary of the trade curve:

* **~6.7 × 10⁴ physical qubits, one run** — if anyone finds a way to produce
  prime-field elliptic-curve relations with known logarithms. 15× current hardware.
* **~2.7 × 10⁷ physical qubits, one run** — as things stand. 6,000× current hardware.
* **Nothing sound in between**, and the smallest piece the problem can be cut into
  is ~9 × 10⁴ qubits.

## 6. Files added in this pass

`instance.py` (verified handle on the core), `weakness.py`, `subsetsum.py`,
`optimality.py`, `embed.py`, `remeasure.py`, `window256_seq.json`,
plus the sequential-counter one-hot in `ladder.py` / `resources.py`.

```bash
python3 weakness.py       # every standard curve weakness, checked
python3 optimality.py     # the cost law and its floor
python3 subsetsum.py      # the encoding this problem wants to be
python3 embed.py          # physical qubits, both worlds
```
