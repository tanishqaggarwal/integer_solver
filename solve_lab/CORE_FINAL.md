# The atom-level obstruction — why `EQUATIONS.txt` (new instance) resists the wire escape

Best verified: **39,022 / 39,033** (`best_agentA_39022.json`). This document records the
**deepest** characterization reached: the residual obstruction, read at the *atom* level with
the exact method that solved the previous instance, and the precise
reason that method is blocked here.

## 1. The obstruction is 3 atoms (wire=p / agent A branch)

Evaluating all 46,298 shared atoms at `best_agentA_39022.json`, **exactly three are nonzero**:

| atom | form | role |
|------|------|------|
| 20862 | `7376877·x_642 + x_2099 − x_7068` | gap **G1** |
| 20864 | `x_4432 − x_19964 − x_28730` | gap **G2** |
| 42669 | deg‑2 `= 24·G1 − 29·G2 + (products that are 0 here)` | verifier tied to G1,G2 |

These three atoms are the entire content of the 11 failing equations
`[2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125]`. (Previous instance had 4
nonzero atoms — same shape.)

## 2. Both gaps are product slacks (the wire-escape structure is present)

- **G1**: the `x_642` term is a multiply gate — `x_642 = x_28599 · x_17325`.
- **G2**: the `x_28730` term is a multiply gate — `x_28730 = x_17499 · x_9413`.

And the escape partners are present exactly as in the previous instance:

- `x_17325`: **free, in only 2 atoms** (its slack def + verifier 42669) — a *rare partner*.
- `x_9413` : **free, in only 2 atoms** — a *rare partner*.

So the raw structure to absorb each gap exists: set `x_642` (via `x_17325`) and `x_28730`
(via `x_9413`) to cancel the gaps, leaving the heavily‑used variables untouched — the exact
move of the wire escape.

## 3. Why it is blocked: the wire is force‑locked to p (the setter's hardening)

The other factor of each slack is a **wire member**:
`x_28599` (81 atoms) and `x_17499` (75 atoms) are both in the **220‑variable identity wire**
(the union‑find class of `x_26064`). In the previous instance the wire root was a **free**
variable `x_15`, so the whole wire could be moved to `sign·1` and the slacks became
fine‑grained (`x_642 = 1·x_17325`), absorbing any gap. **Here the wire root is pinned:**

> atom **37110 = `p − x_26064`** — a pure degree‑1 atom, **no product, no slack** — forces
> `x_26064 = p` exactly, hence every wire member `= sign·p`.

Consequently the slacks are quantised to multiples of p:
`x_642 = p·x_17325`, `x_28730 = p·x_9413`. But the gaps are **sub‑p**:

- G2 needs `x_28730 = x_4432 − x_19964 ≈ 0.47·p` (not a multiple of p) → no integer `x_9413`.
- G1 needs `7376877·x_642 = x_7068 − x_2099` (sub‑p) → no integer `x_17325`.

A granularity‑p slack cannot absorb a sub‑p gap in ℤ. **That is the whole wall.**

### The wire=1 branch does not help
Setting the wire to `sign·1` makes the slacks fine‑grained but lights up **atom 37110**
(`p − 1 ≠ 0`) plus codeword atom **45828**, and 37110 has no slack — it appears in all 13
"unpacking" equations. Agents B and E independently proved this branch integer‑infeasible;
the atom view shows why: 37110 is a slack‑free lock the setter added specifically to defeat
the wire escape.

## 4. Where the real witness lives, and why local search can't reach it

At the true witness, `x_28730 = p·x_9413` forces **`x_4432 ≡ x_19964 (mod p)`** (and
`x_7068 ≡ x_2099 mod 7376877·p`). Agent A's solution sits on the **wrong residue**
(`x_4432 − x_19964 = 0.47·p`). Moving onto the correct residue while keeping the other 16
"ripple" verifier equations satisfied is a global re‑solve:

- The 11 fails are *exactly linear* in the two free leaves `x_4432, x_7068` (rank‑2), and
  those leaves are **independent of the core/loads** — so the conflict is pure ripple.
- Pinning `x_4432,x_7068` to the residue‑correct values and re‑solving the induced **16**
  ripple equations is **first‑order infeasible**: the constrained mod‑p Jacobian stays
  inconsistent as the compensator pool grows (35→67→91 knobs, 35 contradictions fixed), and
  the coupling closure is the whole system (6,215 free inputs / 25,895 equations).

Four independent agents (A global core solve; B whole‑system Newton; C null‑space quadratic;
E product‑slack census) and two further independent methods here all converge on the same
point. The residual is a genuinely hardened trapdoor: crossing it needs the setter's 256‑bit
message (a lattice/CVP problem on the message class) or the secret — not a local repair.

## 5. Tooling added this session
`atom_obstruction.py` (nonzero‑atom finder), `atom_read.py`, `atom_slack.py`,
`wire1_atoms.py`/`wire1_read.py` (wire→1 atom analysis), `atom_eqmap.py`;
`heal_harness.py` (forward‑reconstruct + Jacobian on the wire=p branch),
`an_resid/an_coef/an_diff/an_trace/an_fanout/an_feed/an_closure`, `heal16_solve/gf`,
`heal_grow.py`. Best solution unchanged and re‑verified: **39,022/39,033**.
