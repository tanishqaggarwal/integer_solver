# Squeezing ONE 256-bit modular multiplication (measured, real p = 2^256−2^32−977)

Every encoding passed the faithfulness tests in `verify.py` / `demo_win2.py` before
its size was recorded. Reported triple: **(logical qubits, max clique K, |J| bits)**;
physical ≈ logical × max(1, K/6) on Pegasus, so a K ≤ 6 encoding has no embedding
overhead and the whole problem reduces to minimising column entries.

## 0. The published baseline used the wrong modulus
`resources.py:marginal_window` builds with a **random** odd p (popcount ~125); the
real p has popcount **250**, and the reduction quotient's column load is proportional
to popcount(p). Re-measured honestly a `w=4` wallace window is 900,359 (not 705,907).
All numbers below use the real p; every "baseline" row is `qubo.py` re-measured.

## 1. Best measured general 256×256 modular multiply
Karatsuba depth 4 (leaf ≈ 24) · quotient with p in **non-adjacent form** · 3:2
(wallace) carries → **(99,298 logical, K=5, |J|=2⁶)**, 394,406 couplers, avg deg 7.9.
**99,298 physical ≈ 22.6× a 4,400-qubit Advantage2.** As a squaring: **62,692**.
Baseline was 324,853 (wallace) / 3,330,206 physical (binary) → **3.27× smaller**.

## 2. Pseudo-Mersenne reduction is the biggest single win, and free
`qubo.py` binary-expands the coefficient p in `A·B − C − p·q = 0`; popcount 250 makes
the 257-bit quotient word cost ~64,000 column entries — as much as the entire
partial-product matrix. Writing p in NAF (`2^256−2^32−2^10+2^6−2^4−1`, **6 signed
terms**) collapses it: **−38% logical, −3 bits |J|** (wallace). `mmqb.naf_split` does
this for any constant, so it is not specialised to secp256k1 — a pseudo-Mersenne
modulus simply has a short NAF. (Explicit X_lo + c·X_hi folding is 0.7% worse.)

## 3. Karatsuba/Toom crossover (measured, wallace, best leaf/width)
Karatsuba crosses over at s≈24; Toom-3 at s≈50 and never overtakes Karatsuba. At
s=256: schoolbook 199,931 → Karatsuba 99,298 (0.497×) → Toom-3 107,641. Optimal is
**depth 4, leaf ≈ 18–24**; deeper cuts ANDs but per-node word overhead grows n(3/2)^d.

## 4. Squaring is 50% cheaper and Karatsuba compounds it
AND-cache hit rate is exactly 50% for squarings (transpose symmetry), 0% for general
products. Every squaring in the ladder was already encoded as one; `lam·d`, `lam·m`
share no monomials so nothing more can be shared.

## 5. Carry disciplines
wallace (3:2 tree, one term/column): K=5, |J|=2⁶ — **at the clique floor**. Dadda with
a final ripple does not pay (saves 0.5% vars, K 5→8, +33% physical). Unary/thermometer
carries are strictly dominated (K=781). binary is smaller logically (69,959 unchunked)
but K=103–280, |J|=2¹⁹ → 1.5–4.3M physical.

## 6. Degree reduction: AND ancillas beat Ishikawa 13–100×
Two distinct products in a column share no variable, so every cross term is a distinct
quartic; Ishikawa costs Σ_c C(h_c,2) (quadratic in column height) vs linear for one AND
ancilla. Radix-4 Booth rejected on the same argument (identical AND count + sign glue).

## 7. Classical presolve removes essentially nothing — and that is provable
Persistency ceiling (a preprocessor can only fix variables constant across ALL ground
states), measured by enumeration: ≤3.5% and falling with size, all pure range slack.
Bound/unit propagation to fixpoint on the real 256-bit modmul: **fixes 0 to 27** of
~90–200k variables. Reason: in the ladder each product is another unknown word, so with
A,B free every partial product takes both values over the solution set. Nothing to find.

## 8. The full ladder
Squeezed comb window (w=8): 361,716 (was 1,023,479). **Full 32-window ladder:
3.28×10⁷ → 1.16×10⁷ logical = physical, K=5, |J|=2⁹** — 2.8× better than the published
number, 2,631× an Advantage2, still 4–5 bits over on coupler precision.

## 9. How far the best modmul is from 4,400 qubits — a hard floor
**Not close, and not in a way more effort fixes.** Every partial product must be
linearised before entering an adder, so the AND count is a floor independent of every
carry/reduction/presolve trick:

| | ANDs | vs 4,400 |
|---|---|---|
| schoolbook | 65,536 | 14.9× |
| Karatsuba depth 4 (overall optimum) | 22,711 | 5.2× |
| Karatsuba depth 6 (min AND) | 17,777 | 4.0× |
| recursed to single bits (unattainable, ignores +1 growth/level) | 6,561 | 1.5× |

**Even the unattainable limit is 1.5× the whole machine; the real minimum is 4×, in
partial-product ancillas alone, before any carry/word/reduction bit.** No faithful QUBO
of one 256-bit modular multiply fits 4,400 qubits, and the obstruction is not the carry
discipline, reduction, presolve or degree reduction — all near their floors — it is that
256×256 bit-multiplication has more partial products than the machine has qubits.
Fitting a modmul on Advantage2 needs operand width **s ≈ 40–48** (measured: s=48 → 6,347,
s=32 → 3,019). That is a different (smaller) problem — which is exactly the synthetic
regime this repo now studies.

## 10. Faithfulness
verify.py L0 (full ground-state enumeration, E=0 ⟺ solution), L1 (exhaustive over
inputs for p=13,29,61,127,251,1021,8191), demo_win2.py (whole-ladder, every candidate
scalar): all FAITHFUL, 0 failures. Files: mmqb.py, mm.py, ladder2.py, verify.py,
demo_win2.py, measure.py, crossover.py, window.py, presolve.py, degree.py.
