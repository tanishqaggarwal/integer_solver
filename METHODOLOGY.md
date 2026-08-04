# Methodology — solving the `EQUATIONS.txt` integer trapdoor

This document records how the pure‑integer feasibility problem in `EQUATIONS.txt` was
solved: the structure that was reverse‑engineered, the reasoning path (including the dead
ends), the exact construction that produced the witness, and how it was verified. It is
written to be read top‑to‑bottom by someone who has never seen the problem.

- **Result:** a complete integer assignment for all 38,748 unknowns satisfying **all
  39,031 equations exactly in ℤ**. Deliverable: `SOLUTION.json`.
- **Verifier:** `solve_lab/checker.py` (exact big‑integer evaluation, no floating point).
- **One‑line regenerator:** `solve_lab/build_solution.py`.
- **Constraint honored throughout the final solve:** no SAT/SMT solver — the witness came
  from custom structural analysis, not z3/cvc5.

---

## 1. The problem

`EQUATIONS.txt` is 39,031 lines, each of the form `LHS = 0`, where `LHS` is a fully
parenthesized polynomial over the integers in variables `x_0 … x_38747`. The only
operators are `+`, `−`, `*` (no division, no exponent). Verified structural facts:

- 39,031 equations, all `= 0`; exactly one `=` per line.
- Variable indices span exactly `0 … 38747` (38,748 unknowns); none out of range.
- Task: find integers for every `x_i` so that every equation evaluates to exactly 0 in ℤ.
  Floating point is inadmissible — the numbers reach hundreds of digits.

It is not a random system. It is an **obfuscated arithmetic circuit** — a trapdoor:
a solution exists by construction, but the construction is hidden so that naïve search
cannot find it.

---

## 2. The reverse‑engineered architecture

### 2.1 Shared atoms
The 39,031 equations are built from **46,275 shared sub‑expressions** ("atoms" — gate
residuals). Each atom is meant to vanish at the intended solution. By degree:

| atom degree | count  | role                                             |
|-------------|--------|--------------------------------------------------|
| 1 (linear)  | 20,090 | wiring / identity / negation gates               |
| 2           | 25,468 | multiply gates `x_out − x_a·x_b`, load/div gates |
| 4           | 717    | verifier checks — **perfect squares** `Q²`       |

Each equation is a **random linear combination** of atoms (some equations are a single
perfect‑square atom). Consequence: **if every atom is 0, every equation is 0.** So the
intended solution zeros all atoms; the obfuscation is in making that assignment hard to
find. Atoms are reused ~10× on average, which couples the whole system.

### 2.2 Gate vocabulary
- **Identity / negation** (2‑term linear): `x_a − x_b` or `x_a + x_b`. These rename wires
  — the core obfuscation trick. Union‑find over them collapses aliases (see §4.3).
- **Multiply**: `x_out − x_a·x_b`.
- **Load**: `bit·(x_B − HUGE) − s·x_C` — a boolean `bit` gates a ~290‑bit constant onto a
  wire.
- **Div‑wire**: `x_a·x_out − k·x_b` — defines `x_out = k·x_b / x_a`, with a division‑by‑zero
  trap when `x_a = 0`.
- **Verifier squares**: degree‑4 atoms that are exact squares `Q²` of a degree‑2 form `Q`;
  "check passes" ⟺ `Q = 0` (half the degree).

### 2.3 Degrees of freedom
The active core is one giant residual component (~23.8k variables) driven by ~256 free
control bits plus a few thousand free value‑inputs; the rest is slaved by the gates.

---

## 3. The obstruction — "the twist"

A strong partial assignment (`solve_lab/best/best_partial_39019.json`) satisfies
**39,019 / 39,031** equations. At that point **exactly four atoms are nonzero**:
`{1817, 30378, 40782, 44271}` — the entire obstruction lives here. The 12 failing
equations are precisely the equations whose linear combination includes one of these.

Reading the **raw** equations (not any reformulation) pins the obstruction to two facts:

1. **A conserved invariant gap.** With `G = x_17728 − x_3183`, the partial has
   `G = 63398753350954830538284979531311478224817569395477016427713014637060524103217265241016814 ≠ 0`.
   Across every reformulation this gap relocated but never vanished — it behaves like a
   conserved quantity. This looked like a rigid contradiction "`x_3183` must equal
   `x_17728`, but can't."

2. **The other broken half.** `F0 = 6033033·(x_18274 − x_9770)` is nonzero; the atom
   `1817` demands `x_26977 = F0` but the partial has `x_26977 = 0`.

For a long time the twist read as *rigid*: two big variables forced equal yet provably
unequal. That reading is what made it look like a genuine one‑way trapdoor — and it was
**wrong in a way that turned out to be the whole game** (see §5).

---

## 4. Methodology arc — what was tried, and why it stalled

The campaign ran over many sessions. The honest sequence of approaches:

### 4.1 General‑purpose solving (early sessions)
Structural decomposition + constraint propagation reached **39,013 → 39,019**. GF(2)/mod‑2
and mod‑P propagation confirmed the system is linearly consistent but multiplicatively
dense. SAT/SMT experiments (z3/cvc5, unsat‑core localization, WalkSAT) were run to probe
structure; they localized the nonzero value‑inputs but could not synthesize the witness,
and were set aside per the standing "design your own heuristics" constraint. **No SAT/SMT
was used in the final solve.**

### 4.2 The confluent forward‑evaluator (the productive tool that also misled)
A confluent forward‑evaluator was built from the partial's provenance: orient every gate
as input→output, set the free bits, and propagate. This is excellent for *representing* a
consistent assignment — but it has a fatal property for this problem: it **quantizes both
sides of the twist to coprime units and zeros the slack products**. The witness lives on a
branch the evaluator literally cannot express, so *every* forward‑eval search
(simulated annealing, meet‑in‑the‑middle, greedy, bit‑enumeration, local repair)
plateaued. Much of the mid‑campaign was spent fighting this self‑inflicted blind spot:
chasing control‑bit settings (`x_12779 ≥ 2` to dodge the div‑by‑zero trap), "dirty‑bit"
ripple, and a high‑dimensional co‑activation of thousands of free partner variables.

### 4.3 The identity "wire"
Union‑find over all 2‑term identity/negation atoms collapsed aliases and found **one giant
class of 220 variables** all forced equal (up to sign) to a single free variable `x_15`.
This "wire" is where the escape grounds. Signed union‑find gives each member's sign so
`x = ±V`. Recognizing the wire was the key structural unlock — but under the forward‑eval /
mod‑P framing its significance stayed obscured.

### 4.4 Perfect‑square reduction
The degree‑4 verifier atoms are exact squares `Q²`; `try_sqrt` (`solve_lab/check_square.py`)
recovers the degree‑2 root `Q`. This halves the working degree of every check and is what
later exposed `Q` as the single master constraint (§5.4).

**Net state before the solve:** the structure was fully mapped, but every *search* was
blocked because it was conducted inside the one representation (forward‑eval, mod‑P) that
cannot hold the witness.

---

## 5. The reframing that cracked it

The solve came from **dropping the forward‑eval and mod‑P framings entirely** and working
directly in the raw equation space with exact integers. Five observations, in order:

### 5.1 Only two things are actually broken
Evaluate the atoms at the partial: nonzero set is `{1817, 30378, 40782, 44271}`, which
reduce to two independent facts — `H` and `F` below. Everything else is already 0. So the
task is *not* a global search; it is a **local repair of two residuals** that must not
disturb the 39,019 satisfied equations.

### 5.2 The gap is not rigid — it is a product slack
In the raw text, `(x_17728)−(x_3183)` **never appears alone**: all 16 occurrences are
`((x_17728)−(x_3183)) + x_9982`. The true atom is

```
H = (x_17728 − x_3183) + x_9982
```

So there is no rigid "`x_3183 = x_17728`". The gap is absorbed by `x_9982`, and `x_9982` is
itself a **product**: atom 1818 says `x_9982 = x_12518·x_9897`. Likewise the other half:

```
F = 6033033·(x_18274 − x_9770) − x_26977,   with   x_26977 = x_20510·x_31302  (atom 1816)
```

The "rigid twist" was an artifact of a reformulation that split `H` into two pieces. In
reality both obstructions are **product‑slack activations**:

```
fix H  ⟺  x_9982  = −G   via  x_12518·x_9897  = −G
fix F  ⟺  x_26977 = F0   via  x_20510·x_31302 = F0
```

### 5.3 The multiplier factors live in the quiet wire
`x_12518` and `x_20510` are heavily‑used hub variables (271 and 237 equations) — changing
them naïvely would wreck the system. But signed union‑find (§4.3) places **both of them,
and `x_15`, in the same 220‑variable wire class.** Two exact, integer (not mod‑P) checks
seal it:

- The **signed coefficient sum of the whole wire is 0 in every atom's linear part.**
- Substituting `x_member = sign·V` and holding all non‑wire variables at the partial, the
  coefficient of `V¹` **and** of `V²⁺` is **exactly zero in all 5,233 atoms the wire
  touches.**

So the wire is a **genuinely free parameter**: moving the entire wire to `sign·V` changes
*nothing*, for any `V`. (Its hub members appear in big atoms only in canceling
`+c·x_member − c·alias` pairs, and in product terms only against variables that are 0 at
the partial.)

### 5.4 The rare partners carry the escape
`V` alone is inert (§5.3), so activation must come through the two **rare partners**, each
of which appears in only two atoms — its slack definition and one verifier square:

- `x_9897`: atom 1818 (`x_9982 = x_12518·x_9897`) and the verifier square **atom 40782 = Q²**.
- `x_31302`: atom 1816 (`x_26977 = x_20510·x_31302`) and atom 22049.

`Q` (the 38‑term degree‑2 root of atom 40782) is the master verifier constraint. It
contains exactly the slack terms `−42·x_9982 − x_26977 − x_30323`, the twist terms
`39·(x_3183 − x_17728)` and `6033033·(x_18274 − x_9770)`, plus small linear and wire‑bilinear
terms. Setting `x_9982 = −G`, `x_26977 = F0` and moving the wire drives `Q → 0` (verified),
so the verifier square is satisfied too.

### 5.5 The construction (V = 1)
Everything above collapses to a direct algebraic construction — **no search**:

1. Move the whole wire: each of the 220 members → its sign (`x_12518 = x_20510 = +1`).
2. Set the rare partners: `x_9897 = −G`, `x_31302 = F0`.
3. Set the slack outputs: `x_9982 = −G`, `x_26977 = F0`.
4. Everything else stays at `best_partial_39019`.

This satisfies atoms 1818, 1816, the two twist halves `H` and `F`, and the verifier square
40782 — and, as the exact checker confirms, **all 39,031 equations**. Only **224 variables
change** between the partial and the solution: the 220 wire members plus the four
`{x_9897, x_31302, x_9982, x_26977}`. Notably `x_12779` and `x_24026` stay 0 throughout, so
the entire div‑wire / control‑bit / "dirty‑bit" saga of §4.2 was never on the solution
path at all.

---

## 6. Why it was a trapdoor (why the detour happened)

The construction hides the witness behind exactly the representation a solver reaches for
first:

- **The forward‑evaluator cannot express it.** Orienting the gates and propagating
  quantizes the twist sides to coprime units and forces the slack products to 0 — the
  witness is on the `V ≠ 0` branch, off the evaluator's manifold. Any search built on the
  evaluator is searching a space that provably excludes the answer.
- **mod‑P hides the quietness.** The wire's inertness is an *integer* identity (exact
  cancellation of `V¹` and `V²⁺` coefficients). Mod‑P it looked like a coincidence to be
  re‑checked rather than a lever to be used.
- **The reformulation split the key atom.** Splitting `H = (x_17728−x_3183) + x_9982` into
  two atoms made the gap look rigid, hiding that it is an absorbable product slack. (This
  is also why 2 atoms — 44271, 44272 — remain individually nonzero at the verified
  solution: they are the two halves of `H`, always summed in the real equations, and their
  sum is 0.)

The unlock was methodological: **read the raw equations, work in exact ℤ, and repair the
two residuals locally** instead of searching a global reparametrization.

---

## 7. Verification methodology

The witness was checked adversarially, by paths sharing no evaluation code:

1. **Structural integrity.** `SOLUTION.json` has exactly `x_0…x_38747`, all integer, none
   missing/extra. `EQUATIONS.txt`: 39,031 lines, all `= 0`, operators only `+ − *`, max
   index 38747, no division/powers.
2. **Checker discriminates (not a no‑op).** all‑zeros → 11,679 failing; the solution with a
   single variable (`x_9897`) nudged by +1 → 11 failing; the true solution → 0 failing.
3. **Four independent evaluations, all 39,031 / 39,031:** (a) `checker.py`
   (regex→compile→eval); (b) a deterministic re‑run; (c) an independent regex‑substitution
   evaluator; (d) an **AST‑walk evaluator** that parses each equation with Python's `ast`
   and walks the tree by hand using only `+ − *` — no `eval`, no `compile`, no regex
   substitution.

All agree: every equation evaluates to exactly 0 in ℤ.

---

## 8. Reproduce

```bash
cd solve_lab
python3 build_solution.py            # regenerates best/SOLUTION.json from best_partial_39019
python3 checker.py best/SOLUTION.json # -> satisfied 39031/39031 ... RESULT: OK
```

`build_solution.py` is self‑contained: it rebuilds the signed wire by union‑find, computes
`G` and `F0` from the partial, sets the wire to `sign·1` and the four slack variables, and
asserts the four fixed atoms before writing. Its output is byte‑identical to the shipped
`SOLUTION.json`.

---

## 9. File index

| path | contents |
|------|----------|
| `SOLUTION.json` | the witness — integers for all 38,748 unknowns (root + `solve_lab/best/`) |
| `solve_lab/checker.py` | exact big‑integer verifier for `EQUATIONS.txt` |
| `solve_lab/build_solution.py` | deterministic regenerator of the solution |
| `solve_lab/SOLVED.md` | concise statement of the final method |
| `METHODOLOGY.md` (root) | full methodology and reasoning — this file |
| `solve_lab/check_square.py` | perfect‑square root extraction (`try_sqrt`) |
| `solve_lab/best/best_partial_39019.json` | the 39,019/39,031 partial the solve builds on |
| `solve_lab/NOTEBOOK.md`, `STATE.json`, `RESUME.md` | running research log across sessions |

---

*Bottom line:* the last 12 equations did not fall to a bigger search — they fell to a
change of representation. Once the two broken atoms are read from the raw equations as
product slacks, the 220‑variable identity wire is seen to be exactly quiet, and the two
rare partners are set to `−G` and `F0`, the witness is a four‑step algebraic construction,
verified exactly in ℤ.
