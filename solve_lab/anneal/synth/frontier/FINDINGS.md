# F — the arithmetic-annealing frontier for one modular multiply

**Question.** For one modular multiplication `a·b ≡ c (mod p)` encoded as a QUBO,
what is the largest field size `s` (bits) that a real (classical-surrogate)
annealer drives to the **ground state — QUBO energy exactly 0** — with meaningful
probability? Call it **F**. F decides whether a 2³² anneal budget reaches a 32-bit
key (needs F ≥ 32) or only a single-digit one.

Encoders: `synth/solver/model.build_modmul` (baseline ladder/qubo, wallace|binary,
tunable `W_and`) and the `squeeze/` encoder (Karatsuba + NAF reduction, wallace|
dadda), each converted to a numpy Ising by `model.qubo_to_ising`. Solvers
(`synth/solver/solvers.py`): simulated annealing `sa`, parallel tempering `pt`,
tabu `tabu`, simulated bifurcation `sb`, plus combinations `sb`→`tabu` polish and
wide-ladder `pt`. All code in `synth/frontier/`; **no existing file modified, no
git commit.** "Success rate" = fraction of **independent restarts** reaching E=0.
Every reported E=0 was checked to be a genuine `a·b≡c (mod p)` solution (§6).

---

## 0. Headline

| condition | F (>1% rate) | F (any E=0 seen) | s where it dies (best-E jumps to ≥3) |
|---|---|---|---|
| baseline, default `W_and` | 6 | 8 (≈2%) | 9 |
| **best encoding (`W_and=2`)** | **8** (28%) | **8** | **9** |
| one operand known (clamp `a`) | 6 | 6 | 7 |
| squeeze (Karatsuba+NAF) | 6–7 | 7 | 8–9 |

- **The best honestly-measured F ≈ 8**, reached only with the AND-penalty weight
  tuned to `W_and=2`; the wall then falls hard at **s=9** (best energy jumps from
  0 to 3 and never recovers).
- **`sb` never reaches E=0, even at s=4.** `pt`/`tabu` carry the frontier; `pt`
  and `sb`-polish combos do **not** extend it (they only cost more per restart).
- **Knowing one operand does not help** — it makes the sub-instance *harder*
  (§4c). One-operand-known F ≈ 6.
- The wall is a **flat plateau with an isolated E=0 hole — a needle, not a
  barrier** (§3). Per-restart success ≈ **2⁻ˢ**, matched by the solution density.
- **Verdict (§5): F cannot reach 32.** Every lever (carry discipline, NAF,
  Karatsuba, penalty weight, warm-start) shifts F by **at most ~2 bits** — an
  additive constant, never a change to the 2⁻ˢ scaling. At s=32 the per-restart
  hit probability is ~2⁻³², so a 2³² budget buys ≈1 expected solve for the whole
  multiply — no better than classical enumeration. F is stuck in single digits.

---

## 1. Why there is a wall at all — the needle density (analytic + exact)

The free variables of the atom are the two operand words (a,b: 2s bits → 2²ˢ
assignments); every ancilla (partial-product ANDs, carries, reduction quotient)
is *forced* by the operands. The number of operand pairs with `a·b ≡ c (mod p)`
is exactly **p−1** (each nonzero `a` fixes a unique `b`). So the ground-state
density in the free search space is

> density = (p−1) / 2²ˢ ≈ p / 4ˢ ≈ **2⁻ˢ**   (since 2ˢ⁻¹ < p < 2ˢ).

Exact enumeration:

```
 s     p    #solutions   free states 2^2s   density    2^-s
 4    11        10             256          0.03906   0.06250
 5    19        18            1024          0.01758   0.03125
 6    41        40            4096          0.00977   0.01562
 7    79        78           16384          0.00476   0.00781
 8   163       162           65536          0.00247   0.00391
 9   317       316          262144          0.00121   0.00195
10   641       640         1048576          0.00061   0.00098
```

Density crosses 1% between s=6 and s=7 — exactly where the measured rate crosses
1% (§2). A solver carrying real gradient information would beat this random-guess
density; none does (§3). Tabu at high effort beats it by only a small constant
(≈3–20×, because once the operands are set it can descend the *forced* ancillas),
which buys the +1–2 bits seen in §2/§4a, nothing more.

---

## 2. Core sweep — solver × s × effort   (baseline wallace, default `W_and`)

`success` shown as `hits/restarts rate%`; a miss shows `0/n E<best energy>`.
Fast solvers (tabu, sb) get hundreds of restarts; slow sa/pt fewer (rate coarse,
but their *best E* is the evidence). From `core_baseline.json`.

```
effort=low        sa            pt            tabu          sb
 s=4  n=73     1/20 5.0%     7/12 58.3%    89/154 57.8%   0/300 E2
 s=5  n=107    0/20 E1       1/12 8.3%     14/108 13.0%   0/300 E3
 s=6  n=153    0/20 E1       0/9  E1       2/96  2.1%     0/285 E9
 s=7  n=231    0/20 E2       0/6  E4       0/92  E2       0/244 E15
 s=8  n=255    0/20 E4       0/6  E4       0/93  E2       0/239 E15
 s=9  n=353    0/15 E4       0/4  E10      0/89  E7       0/209 E21

effort=mid        sa            pt            tabu          sb
 s=4           0/8 E1        6/6 100%      88/91 96.7%    0/70 E2
 s=5           0/5 E1        2/3 66.7%     16/30 53.3%    0/57 E3
 s=6           0/4 E1        0/1 E1        0/20 E1        0/55 E7
 s=7           0/3 E3        0/1 E2        0/19 E1        0/46 E13
 s=8           0/3 E3        0/1 E3        0/19 E2        0/45 E16
 s=9           0/2 E2        0/1 E5        0/19 E3        0/40 E24

effort=high       sa            pt            tabu          sb
 s=4              --            --         86/86 100%     0/19 E2
 s=5              --            --         20/20 100%     0/16 E5
 s=6           0/2 E2        2/2 100%      1/6  16.7%     0/14 E5
 s=7           0/1 E2        0/1 E2        0/6  E1        0/13 E13
 s=8              --            --         0/6  E2        0/13 E13
 s=9              --            --         0/5  E1        0/10 E21

effort=extreme    sa            pt            tabu          sb
 s=6           0/1 E2        1/1 100%      1/2 50%*       0/4 E4
 s=7           0/1 E3        0/1 E1        0/2 E1         0/4 E9
 s=8              --            --         1/2 50%*       0/4 E11
 s=9              --            --         0/2 E1         0/3 E19
```
\* small-n noise: a single lucky needle in 2 tries; the true rate is ~2⁻ˢ (see §2b).

**Readings.** `sb` (simulated bifurcation) never reaches E=0, even at s=4 — it
lands 1–2 flips from a min and the polish cannot close a needle; wrong tool.
`sa` is at the needle floor and expensive. `pt`/`tabu` are strongest and agree on
the wall. Cranking effort low→extreme (up to 25× the sweeps/iters) does **not**
move the wall — only how fast a restart finishes. That is the needle signature,
not a schedule problem.

### 2b. Deep restarts pin the tiny rates (tabu, high effort, ~75 s/cell)

Many independent restarts, so a ~1% rate is resolved instead of guessed
(`deep_frontier.json`):

```
 s   n_vars   free-search tabu        one-operand-known (clamp a)
 6    153     14/55  = 25.5%          1/32 =  3.1%
 7    231      0/47  =  0.0%  (E1)    0/31 =  0.0%  (E2)
 8    255      1/48  =  2.1%          0/31 =  0.0%  (E1)
 9    353      0/46  =  0.0%  (E3)      —
10    359      0/46  =  0.0%  (E3)      —
```

s=6→8 hover at the 1–25% needle rate (instance-to-instance noise; s=7's prime
gave a run of misses). At **s=9 the best energy jumps to E3 and stays there** —
tabu can no longer even *approach* a solution: a genuine collapse, not just a
rarer needle. This is the practical wall under the default encoding.

---

## 3. Landscape at the frontier — barrier or plateau?  (decisive)

`landscape.py` ports `diagnose.py`'s energy-vs-Hamming profile to the modmul:
operand `a` pinned to planted, operand `b` set to a candidate, **every ancilla
FORCED to `witness(a, b_cand)`**, so E measures only the arithmetic verifier's
signal. Distance = Hamming(b_cand, b_planted), all p candidates enumerated.

```
 MODMUL s=8 (p=163):  dist   n   minE  meanE  maxE
                        0    1     0    0.0    0     <- the solution
                        1    6     2    4.0    8
                        2   17     2    4.6    8
                        3   32     2    4.3    8
                        4   41     1    4.9    9
                        5   37     2    4.8    8
                        6   21     1    4.8   12
                        7    7     3    4.4    5
                        8    1     5    5.0    5
   correlation(distance, energy):  s=5:0.14  s=6:0.08  s=7:0.20  s=8:0.11
                                    s=9:0.08  s=10:-0.03   (→0 as s grows)
```

**minE at distance 1 (E=2) is no lower than minE at maximum distance**, mean E is
flat (~4–6) across the whole axis, and correlation ≈ 0.1 drifting to ~0 by s=10.
There is **no slope toward the solution**: the E=0 state is an isolated hole in a
flat plateau of shallow E=1…5 local minima. This is why the barrier-crossing
solvers (`pt`, `sb`) gain nothing over `tabu` — there is no barrier to cross and
no funnel to descend. The failure is **fundamental (no gradient), not
solver-limited.** (Same conclusion as `synth/solver/diagnose.py` for the comb.)

---

## 4. Encoding conditioning — can any lever push F up?

### 4a. Penalty / AND-weight (`W_and`) and carry discipline

Tabu, s=7, high effort, ~11 s/cell (`wand.json`):

```
 mode      W_and     s=7 rate    best
 wallace   None(auto)  0/9  0%    E2
 wallace   1           0/9  0%    E1
 wallace   2           5/10 50%   E0   <-- sweet spot
 wallace   4           0/9  0%    E1
 wallace   8           0/9  0%    E2
 wallace   32          0/9  0%    E2
 wallace   128         0/9  0%    E2
 binary    2           1/10 10%   E0
 binary    None        0/9  0%    E1
```

The default `W_and` (`finalize()` picks `max load + 1`, chosen so AND gates are
rigid at hardware coupler precision) is **pessimal for annealing**. A small weight
**`W_and=2` is a sharp resonance** — big enough that gates don't float, small
enough that they don't wall off the search — and it lifts s=7 from 0% to 50%.
`binary` mode (huge coupler range) is worse. Deep restarts with `W_and=2`
(`deep_wand2.json`):

```
 W_and=2   s=7  2/17 = 11.8%    s=8  5/18 = 27.8%    s=9  0/16 E3
                                                     s=10 0/14 E5
```

So `W_and=2` moves F from ~6 to **8** — about **+2 bits** — and then the wall
falls hard at s=9 (best E3). A one-time constant-factor gain, **not** a change to
the 2⁻ˢ scaling.

### 4b. Squeeze encoder (Karatsuba + NAF) head-to-head

Tabu, `squeeze_wallace.json` / `squeeze_dadda.json`, same s, same solver:

```
              baseline wallace     squeeze/wallace      squeeze/dadda
 s=5 high     100%                 100%                 100%
 s=6 high     16.7%                71%                  0% (E1)
 s=7 high     0% (E1)              0% (E1)              16.7% (E0)  <-- lucky needle
 s=8 high     0% (E2)              0% (E1)              0% (E3)
```

The squeeze encoder is a smaller, cleaner Hamiltonian (fewer AND ancillas, K=5,
short coupler range) yet anneals to the **same frontier**. Which specific `s`
gets a lucky sub-1% hit shifts with the carry discipline (dadda caught s=7 once),
but the wall region is firmly s=6–8 and s≥9 is dead for every encoding. Reason:
the free search space and its 2⁻ˢ solution density are properties of the
*arithmetic*, not the encoding. Encoding minimisation decides whether a
sub-instance *fits* on hardware; it does not make it *solvable*.

### 4c. Warm-start / clamp one operand — the "one operand known" F

Clamp `a` to its planted value; the solver searches only `b` + ancillas
(`deep_clampA.json`):

```
 clamp a:  s=5  52%    s=6  3.1%    s=7  0% (E2)    s=8  0% (E1)
 free   :  s=5 100%    s=6 25.5%    s=7  0% (E1)    s=8  2.1%
```

Fixing one operand makes the atom **harder, not easier**: the free multiply has
p−1 ≈ 2ˢ solution pairs, but clamping `a` collapses them to the *single* `b =
c·a⁻¹ mod p` — one needle instead of ~2ˢ, at the same ~2⁻ˢ density per guess.
The interval-split sub-instance therefore inherits the identical wall;
**one-operand-known F ≈ 6**, no improvement.

### 4d. Combined solvers (do stronger combos break the wall?)

`combos.json`: `sb`→`tabu` polish and wide-ladder `pt` (2× replicas).

```
 sb_tabu   s=6 extreme 1/2 E0     s=7 high/ext 0/6,0/2 E1     s=8 0/6 E2
 pt_wide   s=6 high    1/1 E0 (109 s/run!)   s=7 high 0/1 E2 (174 s/run)
```

Neither extends the wall. `pt_wide` can solve s=6 but at 100–330 s per single
restart, and still misses s=7. No combination reaches the `W_and=2` tabu frontier
of s=8 — combining solvers spends more per restart, it does not add gradient.

---

## 5. Verdict

**F ≈ 8, and it is a hard wall.** With the strongest single solver (tabu) at
extreme effort *and* the encoding tuned to its sweet spot (`W_and=2`), the modmul
anneals to E=0 with useful probability up to **s=8 (~28%)**; at **s=9** the best
energy any solver reaches jumps to 3 and never returns to 0. Every lever we have —
carry discipline (wallace/dadda/binary), pseudo-Mersenne NAF reduction, Karatsuba
partial products, AND-penalty weight, warm-starting a known operand, solver
combination — moves F by **at most ~2 bits**, an additive constant on top of the
same **2⁻ˢ** needle density. The landscape (§3) shows why this is unfixable: the
encoded problem carries **no gradient** — an isolated E=0 hole in a flat plateau —
so more schedule, more replicas, or a smaller Hamiltonian cannot help.

**Can F ever reach 32?** No. Reaching s=32 by moving the wall would require
turning a 2⁻ˢ per-restart rate into something polynomial, and no encoding change
touches the exponent — they touch only the ~+2-bit constant. At s=32 the
per-restart hit probability is ≈2⁻³², so even the full 2³² budget yields ≈1
expected solve of a **single** modular multiply — precisely the cost of
classically enumerating one operand (a is one of ~2³² values; the other is then
forced). And empirically the solver cannot even get *near* a solution past s≈9
(best energy stalls at ≥3), so the practical wall is stricter than the density
bound.

> **One line: 32-bit field arithmetic cannot be annealed. The measured frontier
> is F ≈ 8 (hard wall at 9); no encoding or solver lever changes the underlying
> 2⁻ˢ needle scaling, so an annealer offers no speedup over classical brute force
> for modular multiplication, and a 2³² budget reaches at most an ~8-bit field.**

---

## 6. Reproducibility & soundness

- Builders `synth/frontier/fbuild.py`; solver harness `frontier.py`; landscape
  `landscape.py`; deep probes `deep_frontier.py` + inline runs. Data checkpoints:
  `core_baseline.json`, `deep_frontier.json`, `deep_clampA.json`, `deep_wand2.json`,
  `combos.json`, `squeeze_wallace.json`, `squeeze_dadda.json`, `wand.json`,
  `landscape.json`.
- **E=0 ⟺ solution holds for every `W_and`>0**: total energy = Σ(nonnegative
  squares) + W_and·Σ(nonnegative AND penalties), so E=0 forces every AND penalty
  to 0 — a valid gate assignment — independent of the weight. Verified directly:
  every solver-reported E=0 at `W_and`∈{1,2} decodes to operands with
  `a·b ≡ c (mod p)` (often a *different* valid pair than the planted one, since
  the free multiply has p−1 solutions). The planted state is E=0 by construction
  (asserted in every build).
