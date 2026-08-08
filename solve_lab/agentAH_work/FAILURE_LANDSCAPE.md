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

* **`M['live']` has exactly 256 members** and `M['dead']` 128, and — checked against a source
  outside T's model — the 256 values of agent X's `exp2sel` (exponent `i ∈ {0..255}` → selector
  variable, in `agentX_work/xdata.json`, derived independently of `full_model.pkl`) are
  **exactly** `M['live']`.  So `S ⊆ M['live']` with `|S| = n` is a weight-`n` key and
  **`|S|` and `w` are the same axis.**  `|S| = 255` and `|S| = 256` are both legal draws.
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
of `|S|`.**

### The atom count and the equation score are related by an exact, checkable incidence

Agent K's second error was counting atoms where the score counts equations.  The two are not
interchangeable, but on this construction they are related by a closed form that can be checked:

> **score = 39,033 − |{ equations containing at least one still-nonzero atom }|**

verified against the scorer on the fleet's own closures, not fitted:

| closure | nonzero atoms | predicted score | `checker.py` |
|---|---|---|---|
| `close_T32h.json`, `close_M32.json` | 3 (`+ (x3178-(x13720*x21170))`) | 39,033 − 28 = **39,005** | **39,005** ✔ |
| `close_T32g.json` | 4 (`+ (x23514-(x6677*x23504))`) | 39,033 − 37 = **38,996** | **38,996** ✔ |
| all 39,018 closures | 2 (the two targets) | 39,033 − 15 = **39,018** | **39,018** ✔ |

and the *failing index sets* agree elementwise, not just in size.  So on this construction the
atom count is a faithful proxy — but only because the incidence was computed and checked.  Every
score in §3 is still the scorer's, never the proxy's.  A score *below* 39,018 means the closure did not finish; a score above it would
mean the selector set solves the instance.

---

## 3.  The landscape

*(regenerated by `ah_table.py`; raw records in `landscape.json`, one `meta_*.json` and one
`close_*.json` per point)*

All points are **fresh processes**, seeds drawn with `random.Random(seed).sample(M["live"], n)`;
no seed is a prefix of another.  `score` is `checker.py`'s, via its own `load_equations` /
`evaluate_all`.  `THE-15` means the failing set is *exactly*
`[4573,7123,7469,9648,11854,16622,17726,21382,25539,28653,29437,31061,32894,32916,34517]`.

| `\|S\|` | seed | variant | score /39033 | # failing | exit reason | wall s | outers | atoms left | footprint |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 101 | verbatim | **39018** | 15 | `CLOSED_TARGETS_ONLY` | 1 | 1 | 2 | **THE-15** |
| 1 | 101 | fast | **39018** | 15 | `CLOSED_TARGETS_ONLY` | 1 | 1 | 2 | **THE-15** |
| 2 | 101 | verbatim | **39018** | 15 | `CLOSED_TARGETS_ONLY` | 232 | 2 | 2 | **THE-15** |
| 2 | 101 | fast | **39018** | 15 | `CLOSED_TARGETS_ONLY` | 11 | 2 | 2 | **THE-15** |
| 4 | 101 | verbatim | **39018** | 15 | `CLOSED_TARGETS_ONLY` | 191 | 2 | 2 | **THE-15** |
| 4 | 101 | fast | **39018** | 15 | `CLOSED_TARGETS_ONLY` | 16 | 2 | 2 | **THE-15** |
| 8 | 101 | verbatim | **39018** | 15 | `CLOSED_TARGETS_ONLY` | 296 | 2 | 2 | **THE-15** |
| 8 | 101 | fast | **39018** | 15 | `CLOSED_TARGETS_ONLY` | 47 | 2 | 2 | **THE-15** |
| 8 | 101 | fast+bounded | **39018** | 15 | `CLOSED_TARGETS_ONLY` | 37 | 2 | 2 | **THE-15** |
| 16 | 101 | fast | **39018** | 15 | `CLOSED_TARGETS_ONLY` | 53 | 2 | 2 | **THE-15** |
| 32 | 101 | fast | **39018** | 15 | `CLOSED_TARGETS_ONLY` | 214 | 4 | 2 | **THE-15** |
| 32 | 101 | fast+bounded | **39018** | 15 | `CLOSED_TARGETS_ONLY` | 153 | 4 | 2 | **THE-15** |
| 128 | 101 | fast | **38960** | 73 | `MEMORY_BLOWUP` | 1934 | 8 | 6 | 73 eqs |
| 255 | 101 | fast | **38975** | 58 | `OUTER_MAX` | 2139 | 16 | 5 | 58 eqs |
| 255 | 101 | fast+bounded | **38940** | 93 | `TIMEOUT` | 3600 | 14 | 8 | 93 eqs |

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

### A second, `|S|`-independent failure mode of the same routine: a memory blowup

At **`|S| = 24`, seed 101** — a *low* `|S|*, well inside the regime the fleet treats as settled —
the run reached **9.7 GB RSS and was killed by the kernel OOM killer**
(`dmesg`: `Out of memory: Killed process 30381 (python3) ... anon-rss:9727464kB`; the only
kernel OOM event on this box, and it named my process).  The log shows it inside the **joint
two-wire pass**, not the single-wire pass — i.e. in code I did not touch,
`t_close2wj.joint_rootsets` / `tv_roots`:

```python
rs = tv_roots(CFo[a], tw, ma, q, min(e, ex[a]), D)
if ma < m:
    rs = set(b for b in range(m) if b % ma in rs)     # O(m) residues, m = q^e up to 1.66e7
...
out.extend(((b, tw) if flip else (tw, b)) for b in cand)   # then O(|cand|) tuples
```

`m` here is a prime power dividing the atom's handle cofactor, and this instance's cofactors
reach 16,595,977.  One such expansion is ~0.5 GB; the surrounding sample loop can hold several.

**This matters twice over.** First, it is a third way for a high-`|S|` probe to end without
closing that has nothing to do with the instance.  Second, agent T's long-running high-`|S|`
jobs use the same code path, so the failure is live for the fleet, not just for me.

Every AH run after 22:14 sets `RLIMIT_AS` (3 GB by default, `AH_MEMCAP_GB`) so the blowup
surfaces as `reason = MEMORY_BLOWUP` in that run's own record instead of an out-of-memory
event that the kernel resolves by killing whichever process on the box happens to be largest.
The two runs already in flight when the guard was written are covered by `ah_rsswatch.sh`,
which kills by PID at 3.2 GB RSS (it identifies processes by PID, never by command line).

---

## 5.  Verdict

**The decision rule was fixed before the runs, and this is the branch it landed in:
"score ceiling flat in `|S|`, same footprint throughout" — with one qualification about coverage
that I state rather than hide.**

### What the data says

1. **Every run that reached the construction's terminal state scored exactly 39,018 with exactly
   the 15-equation footprint** — at `|S| = 1, 2, 4, 8, 16, 32` here (seed 101, independent of the
   seed-7 chain), and at `|S| = 2,3,5,6,7,8,17,32,64,128` in the fleet's own closures, including
   `|S| = 128` on the **independent** seed 59.  There is **no `|S|` at which the ceiling moves.**

2. **Every score below 39,018 has a named non-instance reason**, recorded by the run itself:
   `MEMORY_BLOWUP`, `OUTER_MAX`, `TIMEOUT`.  Not one run ended with the routine saying it had
   nothing left to try.  In the whole campaign the only exits of the form "no addable collateral /
   nothing moves them" are at `|S| = 8, 17, 32` — **low** `|S|`, in the regime everybody treats as
   settled.

3. **Three distinct routine limits, each of which mimics a high-`|S|` constraint, and none of
   which is one:**

   | limit | where it lives | evidence it is the routine |
   |---|---|---|
   | brute-force modular root enumeration, `O(q^e)` per call with `q^e ≤ 1.66e7` | `solve927g.rootset_pp` | 226 of 232 s at **`|S| = 2`**; replaced by an exact equivalent → **byte-identical output**, 12× faster |
   | materialising a full residue ring / an unbounded `out` list | `t_close2wj.tv_roots`, `joint_rootsets` | 9.7 GB at **`|S| = 24`**, 3.2 GB at `|S| = 96`, `MEMORY_BLOWUP` at `|S| = 128`. `|S| = 24` is *low* |
   | `outer_max = 16` | `t_close2wj.close`'s default | `|S| = 255` was at 5 nonzero atoms and **still strictly decreasing** when iteration 16 cut it off |

   The `|S| = 255` row is the sharpest of the three: under the default cap it reports
   `OUTER_MAX` at 38,975, which read naively is "the ceiling degrades above some `B`".  Raise the
   cap and it keeps descending (155 → 8 nonzero atoms).  **That is what a construction artefact
   looks like, and it is exactly the shape §8 predicts.**

4. **The stall that `UPPER_BOUND_MAP.md` §S5 rests on is resolved and was never evidence.** The
   `|S| = 128` probe it cites ended `NO JOINT ROOT mod 116507 (sampled)` — 400 draws out of
   116,507, i.e. a statement about the sampler.  Agent T's own T36 transposition fix closed it, at
   **two independent seeds**, to 39,018 with the identical footprint.

5. **The footprint invariant is a `grep` result** (§2): those 15 equations are exactly the ones
   containing the two residual target atoms, and the ceiling 39,018 = 39,033 − 15 is what the
   construction produces for *any* selector set that is not the answer.  It cannot vary with `|S|`
   short of solving the instance, so **"the footprint is the same at every `|S|`" was never capable
   of carrying information about `w`.**

### Against the existing closures and the bracket

Any claimed bound `w ≤ B` from this direction must be `B ≥ 128`, because `|S| = 128` closes at two
independent seeds; `B ≥ 128` is **worthless** by `UPPER_BOUND_MAP.md` §8's own table (only
`B ≲ 56` beats rho, only `B ≲ 24` is actionable).  So even the most favourable reading of the
high-`|S|` runs cannot produce a bound with cash value.  Nothing here touches `10 ≤ w ≤ 246`.

### What I did NOT establish

* **No point at `96 ≤ |S| ≤ 255` reached the terminal state *in my runs* inside the budgets used**
  (`|S| = 96` blew memory, `128` blew memory then ran out of clock, `192/224/240/248/252` were
  still queued, `255` timed out at 8 atoms).  The fleet's `|S| = 128` closures carry that regime,
  not mine.  **I am not claiming the routine closes at `|S| = 255`; I am claiming its failure to
  is a clock, a memory cap and an iteration cap, each measured.**
* Seed multiplicity is 1–2 per `|S|` at the values I reached, not the several I intended; the box
  is at load 16–25 on 4 cores and each high-`|S|` run costs 30–60 min.
* Runs still in flight when this was written: worker A (`|S| = 128,160,192,…`) and worker B
  (`|S| = 252,248,240,224,…`), PIDs in `pids.txt`, output in `meta_*_g.json`.  Re-run
  `python3 ah_table.py` to fold them in; the table regenerates from scratch.

### The one-line answer

> **Measured against `checker.py`, the closure ceiling is 39,018 at every `|S|` where the
> construction terminates, with the identical failing set, and that failing set is a syntactic
> consequence of where it terminates.  Every apparent degradation at high `|S|` that this campaign
> has recorded is one of three measured limits of the closure routine — a `O(1.7e7)` root
> enumeration, an unbounded residue-set materialisation that has OOM-killed two processes on this
> box, and `outer_max = 16`.  §8 has no visible effect in the scorer.**

### One thing the fleet should act on regardless of §8

`t_close2wj.joint_rootsets` / `tv_roots` will allocate the full residue ring `q^e` (up to
16,595,977 here) and then an `out` list of `m × |cand|` tuples.  Two kernel OOM kills happened on
this box tonight: mine at 9.7 GB (22:11, PID 30381) and one at **14.6 GB** (23:54, PID 22898) that
was **not** mine — every AH run after 22:14 carries a 3 GB `RLIMIT_AS`, so it cannot reach 14.6 GB.
Any agent running that code path unguarded can take down another agent's job.  A 3 GB `RLIMIT_AS`
costs nothing and converts it into a reportable `MemoryError`.
