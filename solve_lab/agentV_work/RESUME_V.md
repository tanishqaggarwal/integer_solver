# RESUME_V — agent V.

**Two tasks.** §0–§5 are the first: the multi-wire joint solve at `|S| = 17`, which agent T closed
while I was working on it (my contribution there is the independent checker verification of T's
artifacts, an independent reproduction of the residue, and the component/cost structure nobody had
computed). §6 onward is the second, assigned after that redirect: **the structural census of the
handle-less atom population**, which is where the integer side now lives.


**Baseline re-verified by me from cold, first thing:**
`python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
→ `satisfied 39026/39033 (7 failing)`, failing `[12231,12270,12350,14584,18673,22044,29125]`.

Constraints I work under: write only inside `agentV_work/`; other agent directories read-only with
`PYTHONDONTWRITEBYTECODE=1`; no git; long jobs detached with the PID recorded in a file and tested
with `kill -0 $(cat …pid)`, never by command-line matching.

---

## 0. THE SITUATION I WALKED INTO — AND THE FIRST THING I FOUND

My brief said agent T was running a 15-pair two-wire shift on the `|S| = 8` residue and that the
`|S| = 17` joint solve was unbuilt. **That was already stale when I read it.** T had, in the ~20
minutes before I started, gone considerably further: `agentT_work/t_close2w.py` (two-wire fallback
inside the closure loop) and `agentT_work/t_joint17.py` (a three-atom joint solve on one wire pair)
had both run, the latter finishing at 14:07 UTC.

So the first useful thing I could do was not to rebuild what existed but to **verify it**, which
nobody had — `t_joint17.py` prints "dumped close_T17j.json — run checker.py on it" and stops there.
Rule 5 of the ledger exists precisely for this.

### V0 — INDEPENDENT CHECKER VERIFICATION OF T's ARTIFACTS  ✅

| artifact | what it is | `checker.py` | failing set |
|---|---|---|---|
| `close_T2ctl.json` | `\|S\|=2` control | **39,018 / 39,033** | the 15 |
| `close_T7.json` | `\|S\|=7` single-wire | **39,018 / 39,033** | the 15 |
| `close_T8w.json` | `\|S\|=8`, **two-wire** | **39,018 / 39,033** | the 15 |
| `close_T17w.json` | `\|S\|=17`, two-wire pass **stalled** | 39,003 / 39,033 | 30 |
| **`close_T17j.json`** | **`\|S\|=17`, joint solve** | **39,018 / 39,033** | **the 15** |

"the 15" = `[4573, 7123, 7469, 9648, 11854, 16622, 17726, 21382, 25539, 28653, 29437, 31061,
32894, 32916, 34517]` — **byte-identical across all four closing runs**, and identical to the
`|S| = 2` control T established itself.

> **`|S| = 8` and `|S| = 17` both CLOSE over ℤ once multi-wire shifts are allowed, and I verified
> that with `checker.py` rather than with anyone's engine.**  All 927 `c > 1` conditions are
> discharged; the only nonzero atoms left are the two target congruences, which are nonzero because
> the ON-set does not fold to the target — a mod-p fact, not an integer-lift fact.

**This retires check-in 84's headline.** "Closure is a small-|S| phenomenon, the boundary lies
between 5 and 8" was a statement about a **single-wire** solver. It is false as a statement about
the instance. Check-in 86's follow-up ("both outcomes are terminal for the line") resolved in the
direction that moves the boundary, not the one that hardens the negative.

---

## 1. WHAT I BUILT (all in `agentV_work/`)

| file | what |
|---|---|
| `v_base.py` | loads L's calibrated engine from a **private mirror** (`agentV_work/mirror/L`), so nothing outside my directory is written. Validates on load: 3,681 single-handle atoms of 9,032, 927 `c>1`, 2,300 SHIFT wires, p 256 bits — L's and T's own census numbers. |
| `v_poly.py` | **arity-k** exact Newton fit + root finding mod `q^e`, written from scratch (Hensel lift over Cantor–Zassenhaus for large primes, enumeration below 2·10⁵). |
| `v_polytest.py` | self-test of the above **against brute force** at the primes this instance actually uses. |
| `v_diag.py` | reproduces the single-wire `\|S\|=17` end state and censuses the residue's coupling. |
| `v_comp.py` | the **full** component structure of the shares-a-condition graph, both readings, structural and probed. |
| `v_joint.py` | the **general k-wire joint solver** (k arbitrary), component-driven, with both guards. |

### The self-test result, run before anything trusted it

```
roots_pp vs brute force, e=1 : q = 11,19,43,89,199,463,3449,4787,39703,46273  — all ok
roots_pp vs brute force, e>1 : q^e = 27,49,121,361,1849,32,625               — all ok
nd_fit exactness k=1..4      : exact at 40 random points each; degrees [3],[3,3],[3,3,3],[3,3,3,3]
SELF-TEST PASS (0 failures)
```

A root finder that silently *misses* roots turns a negative result into a wrong one, so it is
tested at the sizes it will be used at before it is used. (A root finder that *invents* roots
cannot hurt: every candidate is verified by direct recomputation downstream.)

---

## 2. MEASUREMENTS

* `E.run` on the full 9,032-atom model: **4.9 ms**. `nzcount` (relift + run): **11 ms**.
* Probing influence for all 927 `c>1` conditions: `Σ|wires_of(a)| = 2,683`, so 5,366 `E.run`
  calls ≈ **26 s**. Cheap enough to do exhaustively rather than sample.
* Single-wire `|S| = 17` run (my reproduction of L's): **182 s** wall on a machine at load 7 on
  4 cores — consistent with L's measured 186 s/configuration.

### PROCESS: the command-line-matching failure recurred, in a fourth dress — and so did a new one

The lab's rule is "record the PID at launch, never identify a process by command-line matching".
**Both halves bit me in one launch.**

1. `pgrep -a -f "v_joint.py"` matched **my own shell** (the `bash -c` whose command line contains
   the string) alongside the real job. Fourth occurrence in this campaign after L's two `pkill`
   and one `pgrep`. The rule works; the tool keeps changing.
2. **Recording the PID is not sufficient either.** My launch line lost a newline, so
   `echo $! > job.pid sleep 2` wrote the *wrapper shell's* PID, and `kill -0 $(cat job.pid)`
   dutifully reported **ALIVE** — of a bash process, not of my solver. A liveness test that passes
   on the wrong process is worse than no test.

> **The rule needs a third clause, and this is the one I would add:** after recording the PID,
> **verify it names the job** — `ps -p $(cat job.pid) -o args=` must print your script. Liveness of
> *a* process is not liveness of *your* process. `kill -0` cannot tell the difference; `ps -o args`
> can.

### A wart in `closeS4`/`t_close2w` that changes the reported count

`solve_group3` rejects a trial shift by restoring `vv[w]` — but the internal `nzcount(vv)` it used
to judge the trial calls `relift(vv)`, which **mutates the handle variables**, and those mutations
are *not* restored. So the state after a rejected trial carries stale handle offsets. `closeS4`'s
`close()` happens to end with `relift(vv)` and absorbs them; **a caller that drives the pass
directly and counts without a final `relift` reads 5 nonzero atoms at `|S| = 17` where the true
figure is 3.** I hit this in my first `v_diag.py` run and fixed it rather than reporting the 5.
Recorded because it is exactly the kind of thing that produces a spurious "the count went up".

---

## 3. THE |S| = 17 SINGLE-WIRE END STATE, REPRODUCED  (`v_diag.py 17 V17b`, 183 s)

```
greedy fixpoint            : global nonzero 10, 10 violated c>1 conditions
outer 0  single-wire pass  : 6 accepted shifts
outer 1  single-wire pass  : 0 accepted  -> STALL
NONZERO ATOMS = 3 of 9032
    ((x18956-x37892)-x32237)                c=1          TARGET CONGRUENCE
    ((x24468-x13682)-(12354891*x34243))     c=12354891   TARGET CONGRUENCE
    ((x10261-x8912)-(13040669*x27539))      c=13040669   <-- THE RESIDUE
```

**The residue is `((x10261-x8912)-(13040669*x27539))`, `c = 13040669 = 19·199·3449`, and it has
exactly TWO influencing wires: `x10261` and `x27156`.** Reproduced independently of T (my own
loader, my own mirror, `nonzero atoms of 9,032` with the final `relift` applied) and agreeing with
`agentT_work/t_close2w_T17w.log` atom for atom.

So `|S| = 17` needs **k = 2**, and the pair is **forced, not chosen** — there is only one.

---

## 4. THE FULL COMPONENT STRUCTURE  (`v_comp.py 17 V17`, 65 s)

L reported `[1,1]`. That was two conditions at one point in one search. Here is the whole graph.

### Structural — a property of the instance, independent of any assignment

Wires from the expression plus the value closure, over **all 927** `c>1` conditions:

```
307 components
  conditions per component : {1:67, 2:89, 3:53, 4:50, 5:2, 6:20, 7:20, 8:3, 9:2, 11:1}
  wires per component      : {1:67, 2:13, 3:4, 4:43, 5:10, 6:88, 7:34, 10:1, 11:3, 12:1,
                              13:8, 14:7, 15:4, 17:6, 18:6, 20:2, 21:4, 24:3, 25:1, 27:1, 33:1}
  max = 11 conditions ; max = 33 wires
```

### Probed — influence tested by direct recomputation, at the `|S|=17` greedy fixpoint

```
523 of 927 conditions have NO influencing wire at all and drop out entirely
the remaining 404 fall into 377 components
  conditions per component : {1:363, 2:7, 3:2, 4:4, 5:1}
  wires per component      : {1:347, 2:14, 3:3, 4:1, 5:1, 6:4, 7:2, 8:1, 9:1, 10:1, 13:2}
  max = 5 conditions ; max = 13 wires
```

> **The probed graph is an order of magnitude thinner than the structural one.** 523 conditions
> cannot be moved by *any* shift — they are satisfied and inert. Of the rest, 363 of 377 components
> are singletons. **A wire appearing in an atom's expression very often has zero derivative there**,
> so a coupling census taken from variable sets alone overstates the problem substantially. This is
> the same failure mode the lab has hit repeatedly (a family delimited by the wrong structural
> predicate) and it is why `v_comp.py` reports both graphs rather than one.

### The six components that carry a violated condition at `|S| = 17`

| conditions | wires | moduli |
|---|---|---|
| 5 | 10 | 9025705, 10696593, 10937191, 13097371, 16347615 |
| 4 | 13 | 5930437, 7942211, 10353929, **15333171** ← T's `\|S\|=8` residue lives here |
| 4 | 13 | 1707229, 2264251, 2389855, **13040669** ← the `\|S\|=17` residue lives here |
| 3 | 8 | 1852017, 2992901, 15194385 |
| 2 | 7 | 6672769, 12354891 |
| 1 | 6 | 3849267 |

**T's joint solve used a group of 3 conditions on a 2-wire subset of the third row.** My census says
that component actually spans 4 conditions and 13 wires — so T solved a *subproblem* of it and the
subproblem sufficed. That is a real result and it is also a lucky one: nothing guaranteed it.

---

## 5. WHERE THE JOINT SOLVE STOPS BEING BOUNDED  (`v_cost.py`)

Cost per prime power is `q^(e(k-1))`, so it is governed by the **largest prime power Q in the
component's moduli**, and the affordable component size is

```
k_max(Q, budget) = 1 + floor( log(budget) / log(Q) )
```

At a budget of **2·10⁷ enumeration points per prime power** (~1 minute of `roots_pp` work):

```
   k = 2   affordable while q^e <= 20,000,000  -> 307 of 307 components (100%)
   k = 3   affordable while q^e <=      4,472  ->  37 of 307 components ( 12%)
   k = 4   affordable while q^e <=        271  ->   4 of 307 components (  1%)
   k = 5   affordable while q^e <=         66  ->   0 of 307 components (  0%)
   k >= 5  affordable for NOTHING
```

* largest prime power anywhere in the 927 moduli: **16,595,977**
* median component's largest prime power: **630,067**
* components whose **entire** wire set is jointly reachable: **80 of 307 (26%)** — and those are
  precisely the 1- and 2-wire ones.

> **The number worth knowing is k = 3.** Two wires is affordable for every component in this
> instance. Three wires is affordable for one component in eight. Four is affordable for four
> components out of 307. **Five is affordable for none, at any budget a person would wait for** —
> at `Q = 630,067` (the median) `k = 5` costs `Q⁴ ≈ 1.6·10²³` points.
>
> **This is a hard ceiling and it does not move with cleverness about the search.** It moves only
> if the univariate root-finding is replaced by something that solves the multivariate system
> directly (Gröbner over ℤ/q^e, resultants), which is the one structural idea this line has not
> tried.
>
> **Scope, stated as the rules require:** knob set = the `SHIFT` wires L's `influences()` admits
> (2,300 of them); configuration = the `Random(7)` `|S| = 17` ON-set at its greedy fixpoint, and
> the structural graph over all 927 which is configuration-independent. Budget stated above.
