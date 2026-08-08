# Convex / submodular structure of the modmul QUBO — does it shrink the hard part?

Measured on the schoolbook modular-multiplication QUBO built by `squeeze/mm.py`
(imported read-only) at operand widths s = 8…256, on one comb window, and
exhaustively on small primes. Convention: in QUBO minimisation a coupler is
**submodular** iff its quadratic coefficient c ≤ 0 (ferromagnetic, min-cut
minimisable) and **supermodular** iff c > 0 (the AND/product core).

**Verdict up front:** the ~62 % submodular fraction is real and stable, but it
is a per-*coupler* property, not a per-*variable* block. The supermodular core
touches ≈100 % of the variables, so there is no convex sub-block to peel off by
min-cut. With the operands free the exact persistency ceiling is **0–1
variables** (pure range slack); the only reduction available comes from *pinning
operand bits* (interval-split branching), which is search, not a convex win.
Convex structure **reorganises the arithmetic; it does not shrink the hard part.**
All code in this directory; nothing committed to git.

---

## 1. Submodular fraction vs size — it *holds* at ~62 %

`sweep.py` (modmul), one comb window via `model.build_comb`.

| s | vars | couplers | submodular | supermodular | sub % | sup % |
|---:|---:|---:|---:|---:|---:|---:|
| 8 | 265 | 969 | 609 | 360 | 62.8 | 37.2 |
| 16 | 1 043 | 4 258 | 2 637 | 1 621 | 61.9 | 38.1 |
| 32 | 3 977 | 16 856 | 10 403 | 6 453 | 61.7 | 38.3 |
| 64 | 15 517 | 67 268 | 41 441 | 25 827 | 61.6 | 38.4 |
| 128 | 60 473 | 264 447 | 163 009 | 101 438 | 61.6 | 38.4 |
| 256 | 200 699 | 864 763 | 532 475 | 332 288 | 61.6 | 38.4 |
| comb w=2, μ=8 | 48 615 | 207 519 | 128 712 | 78 807 | 62.0 | — |

The submodular fraction is **flat** (62.8 % → 61.6 %, a slight *decline*). It
neither grows nor shrinks with size, so the convex part never asymptotically
dominates — the supermodular core scales linearly with the problem (∝ couplers).

**Who is the supermodular core.** The positive couplers are the whole
multiplication network, not just the literal AND gates:

| s=256 supermodular couplers | share |
|---|---:|
| adder × adder (Wallace-tree sum·carry `s·d` pairs, coeff +4) | 50 % |
| AND × word / word × word (partial products `a_i·b_j`, coeff +W) | 45 % |
| AND × AND (two partial products in one column, coeff +2) | 10 % |

So "37 %" is not the AND ancillas alone — it is every place two same-sign terms
meet inside a square: the operand-product literals **and** the summation tree
that adds them up. The genuine product core (word×word `a_i b_j`, ~20 %) and the
partial-product interactions (AND×AND, ~10 %) are irreducibly supermodular.

## 2. Persistency — QPBO / probing cannot beat the exact ceiling, which is ≈0

This is a **satisfaction** problem: every penalty is a perfect square / Rosenberg
AND penalty, so all ground states have E = 0. Consequently only **strong**
persistency (a value constant over *all* ground states) reduces the problem;
weak persistency / classical QPBO autarky (a value occurring in *some* optimum)
is vacuous here — both values of any extendable variable occur in some solution.
Strong-persistency roof duality is a **subset** of the constant-over-all-states
set (a theorem), so the enumerated ceiling upper-bounds every QPBO+probing+
autarky method.

Exact, by full ground-state enumeration (each (a,b,c) triple has exactly one
zero-energy completion — the encoding is a deterministic function of the operand
and product words; cross-checked at p=13 against a DFS enumerator, `run_persist.py`):

| p | s | vars | #ground states | **exact ceiling** | propagation | probing |
|---:|---:|---:|---:|---:|---:|---:|
| 13 | 4 | 87 | 351 | 1 | 1 | 1 |
| 29 | 5 | 122 | 1 215 | 1 | 1 | 1 |
| 61 | 6 | 165 | 4 479 | 1 | 1 | 1 |
| 127 | 7 | 190 | 16 892 | 0 | 0 | 0 |
| 251 | 8 | 259 | 67 596 | 0 | 0 | 0 |

The entire ceiling is **one adder carry bit** (`fa:mm:N:7:c`) — pure range slack,
and it vanishes at larger p. Propagation already attains it, so **probing and
QPBO add nothing** (they are bounded above by the ceiling, which propagation
already meets). This is stronger than the prior `squeeze/presolve.py` §7 result:
that used the witness-replay set, which only *upper-bounds* the true ceiling;
here the ceiling is the exact constant-over-all-states set and it is 0–1.

### Conditioned on k pinned operand bits (the interval-split setting)

Exact conditioned ceiling (small p), mean over the 2^k low-bit patterns of
operand a; `[min–max]` over patterns:

| p (s, vars) | k=0 | k=1 | k=2 | k=4 | k=s (a fully known) |
|---|---:|---:|---:|---:|---:|
| 13 (4, 87) | 1 | 5 [2–8] | 10 [3–22] | 24 [5–70] | — |
| 61 (6, 165) | 1 | 6 [3–10] | 12 [3–24] | 25 [5–52] | 44 [7–117] |
| 251 (8, 259) | 0 | 6 [1–12] | 13 [2–26] | 28 [4–61] | 62 [8–184] |

256-bit real modulus, **sound propagation** (a necessary-consequence closure —
sound at any size, verified against the exact ceiling at p≤251; incomplete, so
these are *lower bounds* on the conditioned ceiling), `run_256.py`:

| # pinned low bits of A | fixed (random pin) | fixed (all-zero pin) |
|---:|---:|---:|
| 0 (unconditional) | 1 | 1 |
| 1 | 267 | 267 |
| 8 | 1 810 | 3 968 |
| 32 | 5 671 | 17 339 |
| 128 | 21 242 | 71 839 |
| 256 (A fully pinned) | 45 454 | 173 170 |

Each pinned operand bit cascades to ≈150–260 fixed variables (random pin), more
when the pinned bits are zero (they null whole partial-product rows). Contrast:

* **A and B both fully pinned** → propagation fixes only 47 % (the *true* ceiling
  is 100 % — every ancilla is then determined; the gap is propagation's
  incompleteness, not a fundamental limit).
* **C (the product) fully pinned, A,B free** → 264 / 200 699 = **0.13 %** fixed
  (matches presolve's 266). Knowing the *product* fixes almost nothing — that is
  the factoring-hard direction — whereas knowing *operand* bits cascades.

So the "qubit reduction as a function of answer bits pinned" is real but it is
exactly the interval-split branch cost: reduction ≈ (150–680) × (operand bits
you commit to). Committing operand bits *is* the search.

## 3. Submodular-core decomposition — no separable min-cut block exists

`decompose.py` splits H = H_sub (c≤0) + H_super (c>0) and asks which variables
live *only* in submodular couplers (those are min-cut-determined given the rest):

| s | vars n | free dim 2s | vars in ≥1 supermodular coupler | pure-submodular vars |
|---:|---:|---:|---:|---:|
| 8 | 265 | 16 | 264 | 1 |
| 32 | 3 977 | 64 | 3 973 | 4 |
| 64 | 15 517 | 128 | 15 515 | 2 |
| 256 | 200 699 | 512 | 200 693 | **6** |

**≈99.99 % of variables sit in at least one supermodular coupler.** Only 1–6
isolated variables are purely submodular. The AND core is *woven through the
whole variable set* — even the carry/adder ancillas appear in positive
`adder×adder` couplers (the `s·d` product inside each full-adder square). There
is therefore **no large ferromagnetic sub-block** to solve by min-cut and remove;
submodularity is a property of individual couplers, not of any variable partition.

An actual max-flow min-cut of H_sub in isolation (s=32) returns energy −26 424,
far below the true ground energy 0: minimising the ferromagnetic part alone
neither yields nor bounds a joint solution, because every operand bit also sits
in a positive AND coupler the relaxation discards.

**Are the submodular variables convexly determined by the supermodular ones?**
Functionally yes — the witness is a *total function* of the 2s operand bits, so
all 200 699 variables are determined by 512 free bits (effective dimension = 2s,
0.26 % of the qubit count). But that determination runs through the *entire*
Hamiltonian, not through a convex/min-cut subproblem, and the determining inputs
are exactly the operand bits that feed the supermodular AND core. The convex
structure does not lower the effective dimension; the effective dimension was
already 2s and the search over it is the hard part.

## 4. Honest bound

The convex/submodular structure **reorganises the easy ancillas; it does not
reduce the hard part.** Evidence, all measured and (on small instances)
exhaustively verified:

1. Submodular fraction is stable at ~62 %, so it is not a shrinking obstruction —
   the ~38 % supermodular core scales linearly and spans ≈100 % of variables.
2. With operands free the exact strong-persistency ceiling is 0–1 variables
   (range slack). Roof duality / QPBO / probing / autarky are all bounded by this
   ceiling and propagation already meets it, so **none of them fixes a single
   variable of the multiplication core**. The prior finding ("with answer bits
   free nothing is fixed") is confirmed exactly, and the convex/submodular view
   does not change it.
3. The only reduction is conditioned on *pinning operand bits*: ≈150–680 vars per
   committed operand bit. That is interval-split branching — paid search, not a
   free convex simplification. Pinning the *product* instead fixes 0.13 %.
4. There is no separable submodular block (1–6 pure-sub vars of 200 699), so
   min-cut removes nothing structural. The effective dimension is the 2s operand
   bits, which was already the case; those bits feed the supermodular core that
   the min-cut relaxation cannot touch.

The hard part is choosing the operand bits, which enter through the supermodular
AND core. Convex structure surrounds that core but never removes it.

---

### Files
- `common.py` — read-only imports of `squeeze/mm.py`,`mmqb.py`,`verify.py`; builders + coupler classifier.
- `sweep.py` -> `sweep.json` — submodular-fraction table and supermodular-core kind breakdown.
- `persist.py`, `run_persist.py` -> `persist_small.json` — exact ceiling, propagation, probing; gated by full ground-state enumeration; p=13 cross-checked against DFS.
- `run_256.py` -> `persist_256.json` — sound propagation on the real 256-bit modmul, unconditional + conditioned on k pinned operand bits.
- `decompose.py` -> `decompose.json` — H_sub/H_super split, pure-submodular variable count, and an actual min-cut of H_sub.

### Soundness gate
Every reported fixed variable is constant over the fully enumerated ground-state
set on small instances (`persist.verify_subset`; the replay enumerator is
cross-checked against a DFS enumerator at p=13). Propagation/probing are
necessary-consequence closures (sound by construction) and were verified against
the exact ceiling at p <= 251. No claimed reduction fails the gate.
