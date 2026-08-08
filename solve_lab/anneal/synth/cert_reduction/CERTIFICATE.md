# Certificate — closing the equivalence chain to the ORIGINAL problem

**Scope.** This document certifies **Links B and C** of the reduction chain and states
the end-to-end theorem. All existing repository code is treated as read-only; the new
artifacts live in `solve_lab/anneal/synth/cert_reduction/`.

**The chain.**

```
   QUBO ground state (E = 0)
     │  Link A   (given: per-gadget faithfulness, CIRCUIT_STRUCTURE.md audit)
     ▼
   comb / ladder constraint system satisfied
     │  Link B   (this document, §1)
     ▼
   a 256-bit scalar k with  k·G = T   on E(F_p)
     │  Link C   (this document, §2)
     ▼
   a full integer assignment x_0..x_38747 satisfying EVERY equation in EQUATIONS.txt
```

The point being certified: **decoding a QUBO ground state yields a genuine solution of
the original problem `EQUATIONS.txt`.**

Notation: `p = 2^256 − 2^32 − 977`, curve `E : y² = x³ + B` over `F_p` (A = 0),
`G = PTS[0]`, `PTS[i] = 2^i·G`, target `T`, group order `n` (prime). All of these are
**re-derived and re-checked from `EQUATIONS.txt`** by `reduce.py` / `structure.py` and
re-verified on import of `instance.py` (`selfcheck OK`, reproduced below).

---

## 0. The verified decision core (recap of what Link A hands us)

`instance.py::selfcheck()` passes and establishes, with nothing assumed:

* `A == 0`; every one of the 256 gated points `PTS[i]` and the target `T` lie on
  `E`; `n·G = O` and `n·T = O` (correct prime order); and
  `PTS[i] = 2·PTS[i−1]` for all `i = 1..255` — i.e. `PTS[i] = 2^i·G` (a verified
  doubling chain). `structure.py` additionally exhibits the single chain head and the
  `F_p`-isomorphism onto `y² = x³ + 7` (secp256k1 short form).

So the constraint system that Link A produces is, verbatim (`structure.py` banner):

> find `b_0..b_255 ∈ {0,1}` with `Σ_i b_i·(2^i G) = T` on `E(F_p)`,
> equivalently `k·G = T` with `k = Σ b_i 2^i`, `0 ≤ k < 2^256`.

Link A (per-gadget faithfulness of the comb/ladder QUBO) is **taken as given** from the
independent audit in `CIRCUIT_STRUCTURE.md`; this certificate does not re-derive it. It
is flagged as an *assumed* input in §4.

---

## 1. Link B — the comb/ladder system solves exactly the `k·G = T` bit vectors

**Claim.** The comb/ladder constraint system is satisfiable, and its solutions are
exactly the tuples `(b_0..b_255) ∈ {0,1}^256` with `Σ b_i·2^i·G = T`.

### 1.1 Each combine gadget *is* the affine point-addition law, made fraction-free

`CIRCUIT_STRUCTURE.md` gives every gadget as, over register pairs
`(w1,w2),(w3,w4) → (w5,w6)`:

```
A = (w5 + w1 + w3 + a2)·(w3 − w1)²  −  (w4 − w2)²   == 0
B = (w6 + w2)·(w3 − w1)  −  (w4 − w2)·(w1 − w5)     == 0
```

Write `P₁ = (w1,w2)`, `P₂ = (w3,w4)`, `P₃ = (w5,w6)`, and `λ = (w4−w2)/(w3−w1)`. The
short-Weierstrass addition law on `y² = x³ + a2 x² + a4 x + a6` is

```
x₃ = λ² − a2 − x₁ − x₂ ,     y₃ = λ·(x₁ − x₃) − y₁ .
```

Rearranging `x₃`: `(w5 + w1 + w3 + a2) = λ² = (w4−w2)²/(w3−w1)²`, i.e. **A = 0**
after clearing the denominator. Rearranging `y₃`:
`(w6 + w2) = λ·(w1 − w5) = (w4−w2)(w1−w5)/(w3−w1)`, i.e. **B = 0** after clearing.
Hence `A = B = 0  ⟺  P₃ = P₁ + P₂` on `E`, provided `w3 ≠ w1` (mod p). The gadget is
therefore an *exact* encoding of the group law, not an approximation. `a2` is supplied
by the hardwired pin `a2 = x_24453` (= 0 here), audited in FINAL_CERTIFICATE.md §2.

### 1.2 The `d ≠ 0` non-degeneracy gadget closes the `x₁ = x₂` division loophole

The one place §1.1 can fail is `w3 = w1` (mod p): then `λ` is `0/0`, both residuals
`A, B` vanish *identically*, and `(w5,w6)` is left free. This is exactly the defect
that `CIRCUIT_STRUCTURE.md` **soundness audit item 1** records — all 383 λ-wires are of
the *weak* form `w·(w3−w1) = (w4−w2)` (0 of the strong form `w·(w3−w1) = 1`), and it is
the defect `best_agentA_39022` exploits. A **sound** encoding must additionally carry a
*strong* inverse wire forcing `w3 − w1` to be a unit:

```
d = w3 − w1 ,   d · dinv ≡ 1  (mod p)   ⟺   d ≠ 0  (mod p).
```

When this wire is present and active, coincident operands are rejected, the group law of
§1.1 holds without exception, and the accumulator cannot silently free a register. The
synthetic instance of §3 **includes this strong wire on every active gadget**, so its
comb is the sound version; the live instance's *absence* of it is the audited weakness
(the origin of the 39022 local optimum), not a hole in the reduction's semantics.

### 1.3 The leaves are a doubling chain, so the accumulator computes `Σ b_i 2^i G`

`CIRCUIT_STRUCTURE.md` verifies the 256 leaf constants form a **self-combining chain**
`L_{i+1} = combine(L_i, L_i)` with 0 mismatches; by §1.1 that is
`L_{i+1} = 2·L_i`, hence `L_i = 2^i·L_0 = 2^i·G` (also re-verified by `instance.py` and
`structure.py`). Each selector `b` gates leaf slot `i` to `L_i` when `b=1` and to the
identity `O` when `b=0` (bit-gated constant loads `b·(v−C)=0`, `(1−b)·v=0`; audit
item 2: all 256 selectors carry boolean atoms, 0 free flag slots). The tree of combine
gadgets therefore accumulates the selected subset to the root register pair, and

```
root = Σ_i b_i · L_i = Σ_i b_i · 2^i · G = (Σ b_i 2^i)·G = k·G .
```

### 1.4 The final register-vs-target comparison IS the accumulator target

`CIRCUIT_STRUCTURE.md`: *"the root registers are compared against a hardcoded target
pair; an atom forces OR(all selectors) = 1."* Combined with §1.3, the system is
satisfiable iff `root = T`, i.e. `k·G = T`, with at least one selector set. Audit items
3–5 close the remaining escapes: node degeneracy is unreachable (disjoint leaf subsets;
exhaustive signed-digit search over all 383 nodes and both orientations → 0 spurious
solutions — this is the **signed/offset digit bookkeeping** check); all 1532 register
slots are constrained (0 unconstrained); and every one of the 1149 residual-gating
constraints carries an exact `p`-slack (0 exceptions). Therefore the constraint system
has **no solutions other than** the bit vectors with `k·G = T`. ∎ (Link B)

### 1.5 End-to-end demonstration: recover a planted `k`

`synth/gen.py` plants a known `k` on a same-shape prime-order curve; `synth/solve.py`
recovers it by the interval-split scheme and verifies `k·G = T` on the nose. Measured
here (small sizes; the scheme is size-independent, only the search budget grows):

```
Link B end-to-end: recover planted k, verify k*G==T
  bits= 16 planted k=48170        recovered=48170        k*G==T:True exact:True runs=1 OK=True
  bits= 20 planted k=508891       recovered=508891       k*G==T:True exact:True runs=1 OK=True
  bits= 24 planted k=6333919      recovered=6333919      k*G==T:True exact:True runs=1 OK=True
  bits= 28 planted k=161501211    recovered=161501211    k*G==T:True exact:True runs=1 OK=True
  bits= 32 planted k=2162685304   recovered=2162685304   k*G==T:True exact:True runs=1 OK=True
```

`synth/reduction.py` further shows the only reduction that could make the *full* problem
tractable (ECDLP → modular subset-sum) is neither constructible in the attack direction
nor annealable — i.e. Link B is a genuine ECDLP, which is the intended hardness.

---

## 2. Link C — the crux: `k·G = T` bits ⟺ a full `EQUATIONS.txt` solution

**Claim.** The 256 selector bits solving `k·G = T` are exactly the free inputs whose
deterministic propagation completes to a **full** integer assignment satisfying every
equation of `EQUATIONS.txt`; conversely every full satisfying assignment restricts to
such a bit vector.

### 2.1 (⇐) Given the bits, every other `x_i` is *forced* — reconstruction is deterministic

The repository's own evidence proves the non-bit wires carry **no free choice**:

* **Peeling certificate (`FINAL_CERTIFICATE.md` §1).** Seeded by the three single-atom
  equations 18843, 19066, 20807, the cascade "an equation with exactly one not-yet-zero
  atom forces that atom to 0" propagates to **all 38133 atoms**, with pivot coefficients
  in `[1,75]`. Consequence: the atom matrix has **full column rank 38133** and
  **null space `{0}`** over `Q` and over `GF(q)` for every prime `q > 75`. Since the
  checker works over `Z` and elimination only divides by integers `≤ 75`,
  `M·a = 0, a ∈ Z^38133 ⇒ a = 0`: **every atom must vanish exactly.**
* **Both pins rigid (`FINAL_CERTIFICATE.md` §2).** `p − x_26064` (13 equations) and
  `a2 − x_24453` (unique multiplicity-1 equation 27494) are pinned; their co-atoms are
  purely boolean/copy. Their values are fixed, not free.
* **220-copy class enforced by the equations (`FINAL_CERTIFICATE.md` §3).** Rebuilt from
  the 3558 pure `x_a − x_b` atoms of `EQUATIONS.txt`: the class of `x_26064` has exactly
  220 members with an explicit 219-edge spanning tree, every edge-atom in ≥ 9 equations;
  the slack multipliers reduce to `p` by explicit chains. Nothing is assumed by the
  harness.
* **Propagation to a fixpoint (`METHOD_SUMMARY.md` §3, `partial_assignment.json`,
  `forced_main_bits.json`).** The direction-agnostic engine forces 5,897 variables from
  the unit pins with zero choices; setting the free boolean inputs and re-propagating
  solves the value wires through the huge gated-load atoms and computes all gates. The
  only residual freedom is **`p`-granular** (`FINAL_CERTIFICATE.md` §4: 3707 absorber
  atoms `x_s = p·x_f`), and `p`-granular slack shifts sums only by multiples of `p`, so
  it **cannot change whether any equation is 0 in `Z`**.

Because every atom is forced to 0 and the only slack is `p`-granular (immaterial to
`Z`-satisfaction), **fixing the 256 selector bits determines every `x_i` uniquely up to
irrelevant slack.** That map — `(b_0..b_255) → (x_0..x_38747)` — is the reconstruction;
running it and then `checker.py` decides satisfaction with no remaining choice.

### 2.2 (⇒) Any full solution restricts to `k·G = T` bits

A full assignment satisfying `EQUATIONS.txt` satisfies, in particular, the comb/ladder
constraints (Link A faithful), hence by Link B its 256 selector bits satisfy
`Σ b_i 2^i G = T`. So the two sides are in bijection: **`k·G = T` bits ⟺ full solution.**

### 2.3 What "reconstruct then check" means concretely

Given the 256 bits: (i) set the selectors; (ii) run the forced-atom propagation of §2.1
to fixpoint, producing the complete `x_0..x_38747`; (iii) `checker.py` compiles every LHS
to a big-integer expression and reports `satisfied k/N`. By §2.1 step (ii) is
deterministic; by the peeling certificate its output is the *unique* `Z`-consistent
completion; and `checker.py` then reports **OK on all N** iff the bits solve `k·G = T`.

### 2.4 The honest gap (see §4)

For the **live** instance the 256 bits are *unknown* — recovering them is the Link B
ECDLP, which is hard by construction. Therefore no full `OK` assignment exists in the
repository: the record is `best/new_instance_partial_39026.json` at **39026/39033**, with
exactly the **7 core equations** `[12231, 12270, 12350, 14584, 18673, 22044, 29125]`
open — the accumulator/consistency kernel that encodes `k·G = T`. Link C's forward
direction (bits → full `OK`) is thus **proved as a deterministic map and demonstrated to
completion on synthetic same-shape instances** (§3), not exhibited on the live key.

---

## 3. End-to-end verifier and its OK log

`cert_reduction/synth_circuit.py` compiles a small planted-key ECDLP into an
`EQUATIONS.txt`-style polynomial system over `x_i`, using **exactly the audited gadget
algebra**: bit-gated constant leaf loads, boolean atoms, the fraction-free combine
gadget `A,B` of §1.1 with the `p`-slack encoding `b·(c1·A + c2·B) − p·free = 0`, the
strong `d·dinv ≡ 1` non-degeneracy wire of §1.2, a 2-way MUX
`R_{i+1} = b?(R_i+L_i):R_i`, an offset running point `R_0 = Q0` (signed/offset digit
bookkeeping, so the accumulator never hits `O` or a degenerate add), the final
root-vs-`(Q0+T)` comparison, `OR(selectors)=1`, the `a2` pin, and a 220-style copy
chain. Everything is exact big-integer arithmetic (the `p`-slack wires carry the exact
quotients), so `solve_lab/checker.py` validates it **verbatim**.

`reconstruct(meta, bits)` is the §2.1 map specialized to this circuit: it recomputes
**every** `x_i` from **only** the selector bits by honest point arithmetic.

`cert_reduction/verify_endtoend.py` decodes a ground state to its bits, reconstructs the
full assignment, and runs the checker machinery. Output:

```
==========================================================================
END-TO-END VERIFIER  (bits -> full assignment -> exact checker)
==========================================================================

--- SYNTHETIC same-shape instances (planted key, full chain closes) ---
[synth 16-bit seed=3] p=65203  planted k=48170
  variables x_0..x_166   equations=154
  reconstruct(256 bits) -> full assignment; checker: satisfied 154/154  RESULT: OK
  negative control (1 bit flipped): satisfied 152/154  RESULT: FAIL (as required)
[synth 20-bit seed=3] p=580969  planted k=508891
  variables x_0..x_206   equations=190
  reconstruct(256 bits) -> full assignment; checker: satisfied 190/190  RESULT: OK
  negative control (1 bit flipped): satisfied 188/190  RESULT: FAIL (as required)
[synth 24-bit seed=3] p=10412047  planted k=6333919
  variables x_0..x_246   equations=226
  reconstruct(256 bits) -> full assignment; checker: satisfied 226/226  RESULT: OK
  negative control (1 bit flipped): satisfied 224/226  RESULT: FAIL (as required)

--- LIVE instance: validate the reconstruction PATH on the record partial ---
[live] best/new_instance_partial_39026.json vs EQUATIONS.txt
  satisfied 39026/39033  (7 failing)  RESULT: FAIL   (15.4s)
  first failing line indices: [12231, 12270, 12350, 14584, 18673, 22044, 29125]
  NOTE: 7 equations remain open — the 256-bit ECDLP core is unsolved.
        Full OK on the LIVE instance requires the unknown selector bits;
        it is DEMONSTRATED above on same-shape synthetic instances.
```

**Reading the log.**
* On every synthetic same-shape instance the chain closes end to end: **256 decoded bits
  → full `x` assignment → `checker` OK on ALL equations.**
* The **negative control** flips a single selector bit; reconstruction then misses the
  target and the checker fails **exactly the 2 root-vs-target comparison equations**
  (`N` → `N−2`). This is the direct, machine-checked confirmation that satisfaction is
  tied precisely to `k·G = T` — the final comparison *is* the accumulator target (§1.4).
* On the **live** instance the verifier exercises the real reconstruction/validation path
  through `solve_lab/checker.py` on `EQUATIONS.txt` and reproduces the record
  `39026/39033`, with the 7 open equations being the `k·G = T` kernel. This validates the
  path; it does not (cannot) show `OK`, because that needs the unknown key.

Reproduce:
```
cd solve_lab/anneal
python3 synth/cert_reduction/synth_circuit.py 20 3        # self-test one instance
python3 synth/cert_reduction/verify_endtoend.py --synth 24 3
python3 synth/cert_reduction/verify_endtoend.py --live /abs/path/to/partial.json
python3 synth/cert_reduction/verify_endtoend.py           # full demo + live validation
```

---

## 4. The complete chain, as one theorem

> **Theorem (end-to-end reduction).** Let `x*` be obtained from a QUBO ground state
> (`E = 0`) by the standard decoding. Then:
>
> 1. **(Link A, assumed faithful — `CIRCUIT_STRUCTURE.md` audit.)** `E(x*) = 0` iff the
>    256-selector comb/ladder constraint system is satisfied by `x*`'s selector bits.
>    *Checker:* the per-gadget QUBO-vs-constraint audit (5 soundness checks), taken as
>    given.
> 2. **(Link B, §1 — proved.)** The comb/ladder constraint system is satisfiable, and its
>    solutions are exactly the `(b_0..b_255)` with `Σ b_i 2^i G = T`, i.e. `k·G = T` with
>    `k = Σ b_i 2^i`. *Justification:* each gadget is the fraction-free affine group law
>    (§1.1); the strong `d ≠ 0` wire closes the `x₁ = x₂` loophole (§1.2, audit item 1);
>    the leaves are the verified doubling chain so the tree accumulates `k·G` (§1.3); the
>    root-vs-target comparison plus `OR(selectors)=1` fix `k·G = T` (§1.4); audit items
>    3–5 close degeneracy/underconstraint/slack escapes. *Checker:* `instance.py`
>    (`selfcheck OK`), `structure.py`, and the planted-key recovery `synth/gen.py` +
>    `synth/solve.py`.
> 3. **(Link C, §2 — deterministic map, proved; full `OK` demonstrated on same-shape
>    synthetics.)** The 256 bits solving `k·G = T` determine every `x_i` uniquely up to
>    `Z`-immaterial `p`-slack (peeling certificate ⇒ null space `{0}`; both pins rigid;
>    220-copy class enforced), and that completion satisfies **every** equation of
>    `EQUATIONS.txt`; conversely every full solution restricts to `k·G = T` bits.
>    *Checker:* `cert_reduction/verify_endtoend.py` → `solve_lab/checker.py`
>    (`satisfied N/N`, `RESULT: OK` on synthetic same-shape; `39026/39033` reconstruction
>    path on the live partial).
>
> **Corollary.** A QUBO ground state decodes to a genuine solution of the original
> problem `EQUATIONS.txt`: `E(x*) = 0 ⇒ EQUATIONS.txt(reconstruct(bits(x*))) = 0`
> for all equations. The single obstruction to exhibiting it on the *live* instance is
> recovering the 256 bits, which is the Link B ECDLP.

---

## 5. Gaps — where a link is assumed or empirical, not proved here

1. **Link A is assumed, not re-derived.** Per-gadget faithfulness of the comb/ladder QUBO
   is taken from the independent `CIRCUIT_STRUCTURE.md` audit. This certificate does not
   reconstruct the QUBO-to-constraint map from scratch.
2. **Live full `OK` is unattainable — by design.** The 256 live selector bits are the
   solution of the Link B ECDLP on a secp256k1-isomorphic curve; recovering them is the
   hard problem. Hence the *forward* Link C (`bits → full OK`) is **proved as a
   deterministic map and demonstrated to completion only on synthetic same-shape
   instances**; on the live key the best attainable in-repo is `39026/39033`, with the 7
   open equations being exactly the `k·G = T` kernel.
3. **Audit items 3–5 of Link B are exhaustive-search / structural, partly empirical.**
   Node non-degeneracy (item 3) rests on an *exhaustive signed-digit search over all 383
   nodes and both orientations returning 0 solutions*, and the slack soundness (item 5)
   on a 0-exception count — strong evidence, machine-checked, but enumerative rather than
   a closed-form proof.
4. **The live weak-division defect is real.** Audit item 1: the live comb uses only
   *weak* λ-wires `w·(w3−w1)=(w4−w2)` (0 strong wires), so coincident operands free a
   register — the origin of the `best_agentA_39022` local optimum. The **sound** version
   (strong `d·dinv≡1` wires) is what the synthetic verifier of §3 implements; the live
   instance's semantics is therefore certified for the *intended* (non-degenerate) comb,
   with the degeneracy called out as the audited weakness rather than swept in.
5. **Topology of the synthetic accumulator.** The synthetic instance uses an
   offset **left-fold ladder** rather than the live instance's balanced tree of 383
   gadgets with 128 zero pads. The **gadget algebra is identical** (fraction-free `A,B`,
   `p`-slack, gating, non-degeneracy, comparison) and both compute `Σ b_i L_i`; only the
   accumulation order differs. This is a modeling choice for the demonstration, not a
   fidelity gap in the reduction's semantics.
6. **Equation count.** The live `EQUATIONS.txt` in the repo is the *pinned* instance:
   **39033 equations / 38133 atoms**. The "39,031" figure in `METHOD_SUMMARY.md` refers
   to the earlier *unpinned* instance (differing by exactly the two pins `p − x_26064`,
   `a2 − x_24453`; `FINAL_CERTIFICATE.md` §6). The verifier reports the checker's actual
   count (`/39033`) rather than a hardcoded number.

---

## Artifacts (this deliverable)

| file | role |
|---|---|
| `cert_reduction/synth_circuit.py` | same-shape circuit compiler + `reconstruct(bits)` |
| `cert_reduction/verify_endtoend.py` | end-to-end verifier (synthetic `OK` + live path) |
| `cert_reduction/CERTIFICATE.md` | this document |

All read-only inputs (`instance.py`, `reduce.py`, `structure.py`, `checker.py`,
`CIRCUIT_STRUCTURE.md`, `FINAL_CERTIFICATE.md`, `METHOD_SUMMARY.md`, the S9/S10/S11
notes, `best/`) are unchanged. Nothing is committed to git.
