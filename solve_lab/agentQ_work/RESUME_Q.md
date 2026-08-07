# Agent Q — the algebraic reduction of the instance

This thread's job was the algebra of `EQUATIONS.txt`. It is finished. Everything below was
re-derived **directly from the instance**; nothing depends on another agent's directory, and the
atom database was rebuilt locally with the shared extractor (`qextract.py`, 32,006 gate atoms).
Where another agent's result is used it is named as theirs.

**One-line result.** `EQUATIONS.txt` is satisfiable **iff** `k·G = T` for a 256-bit scalar `k` whose
binary expansion is the leaf ON-set. The remaining work is exactly one discrete logarithm in a
prime-order group with no exploitable structure.

---

## 1. WHAT IS MEASURED

### 1.1 The group
* Substituting `X = x + K/3` (K = `97553848...838891`, present in the instance as the literal wire
  `x_24453`) removes the universal offset and turns the stage law into the plain chord construction.
* Fitting `Y² = X³ + aX + b` from two leaf pins gives **a = 0**, and **256/256 leaf points lie on
  it** — non-singular.
* The law is **measured, not assumed**: associative **297/297**, commutative **297/297**, and equal
  to the shifted chord law **198/198** on random triples.
* Order by Cornacchia on `4p = L² + 27M²`, verified by exact scalar multiplication on 5 points:
  `N = 115792089237316195423570985008687907852837564279074904382605163141518161494337`,
  a **256-bit prime**. `N ≠ p` (not anomalous); `p^k ≠ 1 mod N` for k ≤ 24 (no small embedding
  degree). Trace `t = p + 1 − N = L`.
* `lam.py` confirms the CM endomorphism: with β a cube root of 1 mod p, `φ(X,Y) = (βX, Y)` **is**
  multiplication by a cube root λ of 1 mod N. It buys **√3 only** and does not move the 2¹²⁸ figure.

### 1.2 The leaves
* Pin atoms have the shape `(x_g)*((x_w) − BIGCONST)`. Scanning for them gives **256 selectors,
  each with exactly 2 pins**. With the correct shift (`X = c + K/3 mod p`) **256/256 lie on the
  cubic**.
* Doubling closes into a **single chain of length 256**, and **256/256 satisfy `L_i = 2^i·G`**,
  with `G` the leaf of selector `x2779`. Nothing is inferred.
* The root target is the instance's only pair of unconditional >60-digit pins:
  `T = (C1 mod p, C2 mod p)`, on the same cubic (the swapped orientation is not). Root pin wires
  are `x_24468` (X) and `x_18956` (Y).

### 1.3 The stage gadgets — 383/383
The stage gadget is the division-free chord law (wires hold raw `u = X − K/3`; the `+K` is `3·(K/3)`):

```
dx = ua−ub    dy = ya−yb
R1 = S*dx² − dy²      S = u3+ua+ub+K        <=>  lambda² = u3+ua+ub+K
R2 = A*dx  − B*dy     A = y3+yb, B = ub−u3  <=>  y3+yb = lambda*(ub−u3)
```

Searching the atom DAG finds **383 gadgets**: 89 leaf-adjacent, 78 mixed, 216 internal. Each was
tested by **Schwartz–Zippel on random curve points** against the *actual* sub-DAG:

> **383/383 verified, including 89/89 leaf-adjacent, all with orientation (+1,+1)** — every stage
> computes the plain sum `P_a + P_b`, no sign flips. None of the 1,532 stage core wires is
> multi-defined, so the test used the real gate relations.

**Census — the 178|78 split, derived from gadget arity.** Counting hard-zero input wires:
89 leaf-adjacent (0 zeros, combine two leaves, consuming **178** leaves), 78 mixed (2 zeros — one
leaf plus a dummy: pass-throughs, consuming the remaining **78**), 191 live internal, 25 dead.

### 1.4 The mux layer — 383/383
Each slot is a three-way weighted mux. Read verbatim off the instance at one slot (leaves `2^0`,
`2^164`; selectors `x_2779`, `x_34715`, both boolean-pinned):

```
cA = a(1−b)    cB = b(1−a)    cC = a·b
Xout = cA*Xa + cB*Xb + cC*u3        Yout = cA*Ya + cB*Yb + cC*y3
live_out = (a+b) − ab = a OR b
```

Evaluated on the real leaf constants, all four quadrants match: `(0,0)`→identity `(0,0)`;
`(1,0)`→leaf `2^0`; `(0,1)`→leaf `2^164`; `(1,1)`→**the sum**. Generalised with both coordinate
muxes required to use the *identical* coefficient wires (which rules out accidental matches):
**383/383, zero unmatched** (40 with boolean-pinned live bits, 343 with internal ones).

The identity value `(0,0)` is not a curve point, but it is only ever **passed through** — it can
never enter a chord, because `cC = ab = 0` whenever a child is dead.

### 1.5 The tree
* every one of the 383 slots emits `OR(s1,s2)` of its own two live bits — **383/383**;
* those ORs give **382 parent←child edges** among 383 slots — exactly a tree;
* **exactly one slot has no parent**, and **all 383 are reachable from it**;
* **all 256/256 leaf selectors appear under that single root**;
* the 766 live-bit slots decompose as **256 leaf selectors + 382 child ORs + 128 hard zeros**
  (the zeros being the dead dummy branches of the pass-through slots). Nothing unaccounted for.

### 1.6 The coordinate hand-off — **mod p**
A slot's mux output `M` is not literally its parent's input `P`; they are related by an **affine
alias** `P = M + (multiple of) Q`, with **all 575 slack wires products of two wires**. Of 766 mux
outputs, **573 alias to a parent slot input and 2 alias to the ROOT PIN**
(`x_24468 = x_13682 + 12354891·x_34243`).

**Agent L supplied what I could not see:** the shared slack factors are the constant **p** (220 such
wires in the instance) — which is why nothing forces them to zero and why my 523/52 split was one
population, not two. Every slack term is `p·(free variable)`, so **slack ≡ 0 (mod p)** and
**`parent_input ≡ mux_out (mod p)` unconditionally**. The hand-off follows the tree of §1.5
**mod p**. Over ℤ the slack is free, and its residue is the 927 `c > 1` divisibility conditions,
which are L's and are open.

### 1.7 All-atoms-zero is forced
Agent F's incidence matrix (39,033 × 39,033, 525,982 nonzeros) has **rank 39,033, `ker(M) = 0`**, by
three independent computations; agent T re-verified the certificate from cold and, separately,
established **faithfulness** of M by exact list equality against `checker.evaluate_all` at 10 points.
So in that decomposition all-atoms-zero is an **equivalence**.

### 1.8 The reduction, and what the searches now say
Chaining §§1.2–1.7: **a satisfying assignment with ON-set S ⟹ `k = Σ_{i∈S} 2^i` satisfies `k·G = T`.**
Every link of that implication is a point identity **mod p**, which is the modulus this direction
needs; the 927 conditions over ℤ sit on the *converse* (existence) direction. The searches are
therefore statements about **the instance**:

| program | family excluded | status |
|---|---|---|
| `dlp_bsgs.py` | k < 2⁴⁴ and N−k < 2⁴⁴ | **excluded** |
| `lowwt.py` | Hamming weight(k) ≤ 6 | **excluded** |
| `window.py` | all ON-bits inside a 34-bit window | **excluded** |
| `smallmul.py` | m·T on the ladder, m ≤ 10⁷ (k = 2^i/m mod N) | **excluded** |
| `lam.py` | k = ±λ^j·2^i ; k = a+bλ, \|a\|,\|b\| < 2²¹ | **excluded** |
| `wt7.py` | Hamming weight(k) ≤ 7 | **33.7% covered, no hit — a partial, not a bound** |

### 1.9 What the 39,026 deliverable actually does
Wire census mod p over all 38,748 wires: `2^72·G`'s x-coordinate on **92** wires, `2^235·G`'s on
**5**, their group sum on **0**, the target C1 on **4**. Its ON-set `{2081, 24601}` is ladder
exponents **{72, 235}**. It **does not fold**: it passes a single leaf up the tree as a chain of
one-live-input pass-throughs, cuts the second leaf after 5 wires, and overwrites the value with the
target near the root, paying 7 broken atoms. That is the entire content of the partial.

### 1.10 Two consequences worth carrying
* **The atom is not the unit of failure.** `ker(M) = 0` forbids *all* equations holding with some
  atom nonzero; it does **not** forbid an atom being nonzero inside an equation that still sums to
  zero. The deliverable is exactly that case: six of its seven failing equations contain only atoms
  occurring in **6–15** equations, yet only 7 break — so most occurrences **cancel**. Compensation
  is already happening in the lab's best assignment, and the gap runs in the **favourable**
  direction. (Routed to M.)
* **The degenerate branch is vacuous, not doubling.** Feed a gadget two equal live points: `dx = dy
  = 0`, so `R1` and `R2` vanish **whatever the output is** — **383/383, even with the output set to
  a random wrong value**. The circuit does not implement doubling. The fold picture needs no two
  equal points ever to meet; children of a slot sum over **disjoint** leaf subsets, so they coincide
  only if `Σ_{S1} 2^i − Σ_{S2} 2^i = ±N`, both sums being `< 2²⁵⁶ < 2N`. That is a checkable
  condition on the particular scalar, not a generic hazard. (Same criterion K and the coordinator
  reached by two other routes.)

---

## 2. WHAT I RETRACTED, AND WHY

**§14(a) — the strong reading of the propagation result.** I measured that setting a selector ON
puts that leaf's coordinate on **no wire at any weight tested** (0/1, 0/2, 0/3, 0/5, 0/7, 0/128;
520/8,583 free inputs solved; 0 contradictions) and reported it as *routing is a constraint, not a
propagation*. That much stands. But I let it read as an absence of determination, and it is not:
a leaf pin is `sel·(w − C) − z`, not `sel·(w − C)`, so the coordinate lands only once `z` is
separately forced to 0. **Routing *is* determined — by a simultaneous system, not by propagation.**
The table measures the weakness of unit propagation. This also partly reconciled agent P, whose
"liveness is determined by the selectors" was closer to right than the contradiction suggested; the
disagreement was about *how*, not *whether*.

**§22 — the atom-forcing gate I opened.** I found 47,198 distinct atom terms against 39,033
equations and argued the incidence matrix has nullity ≥ 8,165, so the bundling does not force atoms
to zero. **That was an artefact of my own parser.** Of my 47,198 terms, **39,032 occur in ≥2
equations — matching F's 39,033 to within 1 — and 8,166 occur in exactly one, accounting for exactly
the 8,165 excess.** My parser splits sub-expressions F treats as single atoms, each split landing in
one equation. I had flagged in §22 that my finer atoms are constrained functions of wires rather
than free coordinates; that caveat was the whole story. I should have run the reconciliation before
reporting the gate. With F's `ker(M) = 0` and T's faithfulness, the gate is closed.

**§9 → §15 → §24 — the six search programs, withdrawn and restored.** Their journey, with the
reason at each step:
1. **Reported with full standing.** Six negative sweeps over structured families of `k`.
2. **Withdrawn (§15), unprompted.** They computed the fold **inside the group model** and never
   checked the circuit agrees — and my §14(a) measurement showed the circuit-side check does not
   close by propagation, in exactly the low-weight regime where `lowwt`/`wt7` live. Their clean-miss
   verdicts were therefore evidence about the **group model**, not the instance.
3. **Held through two check-ins**, including one where the mux law was confirmed at one slot and
   another at 188/383 — neither of which was enough.
4. **Restored (§24)** only once every link of *assignment ⟹ scalar* was measured: all-atoms-zero
   forced, leaves, gadgets, muxes, tree, and the hand-off mod p — with the argument that **mod p is
   precisely the modulus the negative direction needs**, the 927 over ℤ sitting on the converse.

The searches themselves never changed. What changed was what could honestly be claimed from them.

---

## 3. THE TWO RULES THIS THREAD PRODUCED

**Rule 1 — a count derived from one parse is a fact about that parse until reconciled.**
This lab has **five** atom counts: **39,033** (F), **39,277**, **40,727**, **40,885** (I),
**42,267** (A/G/H). Kernel dimension differs across them: 0 in F's, ≥1,852 in I's, ≥3,234 in
A/G/H's. My own 47,198-term parse produced a nullity bound that evaporated on reconciliation, and
my earlier "every atom appears in exactly 1 equation" was an artefact of reading the *deduplicated*
`gates.jsonl`. Agent O hit the same class of error from the other side — `eq8680` has one term in
H's parse and twenty in E's — and caught it by cross-checking rather than by reasoning harder inside
one parse. **Any claim resting on a count must name its decomposition, and must be reconciled
against at least one other before it is load-bearing.**

**Rule 2 — decline to close on structure alone.**
Twice I had a result whose *shape* was exactly what was needed and stopped short of claiming it:
* the coordinate hand-off had precisely the affine-alias form that would close the existence result,
  and I would not claim it because I had measured the tree on the **liveness** side and not the
  coordinate side. L's constant-p result later supplied the missing measurement — mod p.
* the quadrant law held at one slot and then at 188/383, and I would not generalise from an
  isomorphism I had not measured; the remaining 195 turned out to be **my own labelling artefact**
  (`qstages.py` assigned `(u3,y3)` as `sorted(free)[0],[1]`, which is X/Y order at some slots and
  reversed at others), and fixing it *while tightening the test* gave 383/383.

Both times the structure was right. Both times the claim would have been unearned when made.

---

## 4. WHAT IS OPEN, AND WHOSE IT IS

* **The 927 `c > 1` divisibility conditions over ℤ** — the existence direction. **L's**; |S| = 2 is
  closed and |S| = 17 is one step away.
* **The collision criterion** `Σ_{S1} 2^i − Σ_{S2} 2^i = ±N` — needs the particular scalar.
* **Atom compensation as a scoring lever** — **M's**, with the exact scorer.
* **The discrete logarithm itself** — ~2¹²⁸ generic, no exploitable structure found: prime order,
  non-anomalous, no small embedding degree, √3 from the endomorphism. Nothing about the circuit
  reduces it, because the fold is a group homomorphism of the selector vector.

Nothing on the circuit side remains for this thread.

---

## 5. ARTIFACT INDEX (all in `solve_lab/agentQ_work/`)

| file | what it does |
|---|---|
| `qextract.py` | rebuilds the atom database locally from `EQUATIONS.txt` |
| `qgrp.py`, `fastg.py` | group law on the cubic; fast Jacobian + gmpy2 + batch inversion |
| `qpins.py` → `qleaf.json` | all 256 leaf pins straight from the instance |
| `qladder2.py` → `qladder.json` | the single 256-long doubling chain, `L_i = 2^i·G` |
| `qstages.py` → `qstages.json` | 383 chord gadgets, Schwartz–Zippel verified |
| `qdegen2.py` | the degenerate branch is vacuous (383/383) |
| `qmux2.py` | quadrant mux law at 383/383, both coordinates, same coefficients |
| `qquad.py` | the four quadrants numerically at one slot |
| `qlivetree.py` | the liveness tree: 382 edges, one root, 256/256 selectors |
| `qalias.py` | the affine-alias hand-off layer |
| `qsolve.py`, `qrun2.py` | unit propagation over all 47,198 terms; the non-circular routing test |
| `qmult.py` → `qmult.pkl` | true atom/equation multiplicities (non-deduplicated) |
| `dlp_bsgs.py` `lowwt.py` `wt7.py` `window.py` `smallmul.py` `lam.py` | the six scalar searches |

No infeasibility is claimed anywhere in this thread. No "nothing can move X" claim is made. No
named-object framing and no generator forensics were used at any point.
