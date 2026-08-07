# FAILURE_LANDSCAPE — the closure ceiling as a function of `|S|`, measured with `checker.py`

Agent AH.  Working directory `solve_lab/agentAH_work/`.  `PYTHONDONTWRITEBYTECODE=1`.
No git commands were run.  No file outside this directory was created, modified, moved or deleted.

**Verified before anything else:**
`python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json` →
`satisfied 39026/39033  (7 failing)  [12231, 12270, 12350, 14584, 18673, 22044, 29125]`, `RESULT: FAIL`.
My tooling agrees with the fleet's.

---

## 1.  Which construction I used, and what I checked about it

I did not write a construction.  I used the fleet's, unchanged:

```
t_close2wj.close(S, ...)                 agentT_work/t_close2wj.py     (agent T, "AUDIT T33")
  └ t_close2w.two_wire / fit2 / probe2   agentT_work/t_close2w.py      (agent T, "AUDIT T32")
     └ closeS4.solve_group3 / nzcount    agentT_work/mirror/L/closeS4.py   (agent L)
        └ closeS3.roots_c / fitc         agentT_work/mirror/L/closeS3.py
           └ solve927g.rootset_pp        agentT_work/mirror/L/solve927g.py
              └ solve927.fit / peval     agentT_work/mirror/L/solve927.py
                 └ mkassign2 / calib2.assignment / ORIENT / PIN / M['live']
```

`ah_run.py` reproduces `t_close2wj.close` **line for line** and adds only instrumentation
(§2).  Everything it calls — `assignment`, `relift`, `solve_group3`, `joint_pair`,
`handleless_pass`, `forced_exact_pass`, `nzcount` — is imported from `t_close2wj` /
`t_close2w`, not re-implemented.

What I checked about it before trusting it:

* **`M['live']` has exactly 256 members** and `M['dead']` 128.  A selector set `S ⊆ M['live']`
  with `|S| = n` is the construction's stand-in for a weight-`n` key, so `|S|` and `w` are
  the same axis.  `|S| = 255` and `|S| = 256` are both legal draws.
* **The guard is by direct recomputation, not by the fit.**  `nzcount` calls `relift` and a
  full `E.run`; `solve_group3` accepts a shift only if the *global* nonzero-atom count
  strictly decreases; `joint_pair` re-probes every atom in the group with `probe2` before
  accepting, and validates each fitted bivariate against `probe2` at three points.
  This is the discipline agent K's routine lacked (it solved constraints backwards).
* **The one pass allowed to raise the count is `forced_exact_pass`**, and only for a
  handle-less atom whose residual is *linear* on its single admitted wire, so the root is
  unique; it re-checks `E.run(vv)[i] == 0` after applying and freezes the wire.  I left it in.
* **`nzcount` counts ATOMS (9,032), the scorer counts EQUATIONS (39,033).**  These are not
  the same number and I never report one as the other.  Every number in §3 came out of
  `checker.py`'s own `load_equations` / `evaluate_all` (`ah_score.py`, `ah_table.py`), which
  reproduce `checker.py`'s verdict on `new_instance_partial_39026.json` exactly.
* **Per-run process isolation.**  `forced_exact_pass` mutates the module-global `SHIFT`
  (`SHIFT.discard(w)`), and `t_close2w.rnd` is module-global.  Two closures in one process
  would contaminate each other, so every point in §3 is a **fresh process**.
* The historical `|S| = 2` special case (`S = [24601, 2081]` when `seed == 7`) is honoured
  only for seed 7; every other seed draws `random.Random(seed).sample(M['live'], n)`.
  **No seed in §3 is a prefix of another** — that is the specific defect in the existing
  `|S| = 32/64/128` data and it is not reproduced here.

---

## 2.  What the "identical 15-equation footprint" actually is

Every closure the fleet has reported at 39,018 fails exactly

```
[4573, 7123, 7469, 9648, 11854, 16622, 17726, 21382, 25539, 28653, 29437, 31061, 32894, 32916, 34517]
```

and ends with exactly two nonzero atoms,
`((x18956-x37892)-x32237)` and `((x24468-x13682)-(12354891*x34243))`.

**These are the same fact.**  Grepping `EQUATIONS.txt` for the two atoms as literal
subexpressions:

| atom | equations containing it |
|---|---|
| `(((((x_18956)-(x_37892))))-(x_32237))` | 14 |
| `(((((x_24468)-(x_13682))))-((12354891)*(x_34243)))` | 13 |
| **union** | **15 — exactly the footprint above** |

So "the identical 15-equation footprint across ten `|S|` values" carries **no information
beyond "the closure terminated in the same state"**.  It is not an invariant of the instance;
it is the syntactic incidence of two atoms.  Anyone reading it as a mysterious constant is
reading a `grep` result.

And those two atoms are the two the construction *pins*: `assignment()` sets
`v[24468] = T1`, `v[18956] = T2`, the target's coordinates.  The two residual atoms are the
statement **"the ladder output equals `T`"** — i.e. the ECDLP condition itself.  They are
nonzero for every `S` that is not the answer, which is why 39,018 is where this construction
stops.  **39,018 = 39,033 − 15 is the construction's structural ceiling, not a measurement
of `|S|`.**  A score *below* 39,018 means the closure did not finish; a score above it would
mean the selector set solves the instance.

---

## 3.  The landscape

*(regenerated by `ah_table.py`; raw records in `landscape.json`, one `meta_*.json` and one
`close_*.json` per point)*

TABLE-PLACEHOLDER

---

## 4.  Why the high-`|S|` probes looked like a constraint — measured, not guessed

**This is the part the fleet's evidence could not distinguish, and it is a property of the
routine.**

The single-wire pass reaches the instance through `closeS3.roots_c`, which calls
`solve927g.rootset_pp(coeffs, q, e)`:

```python
def rootset_pp(coeffs,q,e):
    m = q**e
    return {t for t in range(m) if peval(coeffs,t,m)==0}
```

— brute force over the **whole** residue ring.  Measured over the instance's own handle
cofactors (`SL[a]/p`):

* **927 atoms carry a cofactor `c > 1`**;
* the largest prime power dividing such a `c` has **median 19,717 and maximum 16,595,977**;
* **19.2 %** of those atoms carry a prime power above 10^6 and **4.5 %** above 10^7.

One call can therefore cost 1.7 × 10^7 big-integer Newton evaluations.  Measured directly:

| point | wall | inside `solve_group3` | `solve_group3` calls | guard (`nzcount`) evals |
|---|---|---|---|---|
| `\|S\| = 2`, seed 101 | 232 s | **226.6 s** | 12 | 5 |
| `\|S\| = 8`, seed 101 | 296 s | **280.6 s** | 26 | 253 |
| `\|S\| = 160`, seed 101 | killed at 469 s | 277 s | 11 | 19 |

At `|S| = 160` the run was at **wire 10 of 439** after 469 s (`evidence_n160_plain_profile.txt`),
i.e. outer iteration 0 alone projects to ~1.1 × 10^4 s.  Nothing about that is the instance:
it is `≈25 s` of modular root enumeration per wire, and the number of wires grows linearly
with `|S|` (26 wires at `|S| = 8`, 262 at 96, 439 at 160, ~600 at 250).

**A construction whose cost is linear in `|S|` with a 25-second constant looks exactly like a
constraint at high `|S|`, and is not one.**

### The fix, and why it is not an algorithm change

`ah_roots.py` replaces `rootset_pp` with the **same set** computed by agent T's own
`t_poly.roots_pp` (root-find mod `q`, then Hensel-lift to `q^e`) — the routine
`t_close2wj` already uses for exactly this purpose in its two-wire pass.  Validation:

* `ah_roots_selftest.py`: **1,984 comparisons against the brute-force original**, on the
  919 distinct `(q,e)` actually occurring in this instance's cofactors (all 491 with
  `q^e ≤ 30000`, four random degree-4 coefficient draws each, plus identically-zero cases).
  **0 mismatches.**  It exercised all four code paths (864 fast, 1,115 brute, 4 `newton_to_mono`
  fallbacks, 1 `ALL`).
* End-to-end: at `|S| = 4` and `|S| = 2` the closure output is **byte-identical**
  (`md5 7da101fc…` and `md5 9f89ffd9…` for plain and fast), at 16 s vs 191 s and 24 s vs 232 s.

Every run in §3 marked `fast` uses it.  Runs marked `plain` do not.

---

## 5.  Verdict

VERDICT-PLACEHOLDER
