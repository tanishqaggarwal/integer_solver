# Agent Y — the complement identity, and the centre framing that subsumes it

**Angle complete.** This file is the consolidated hand-off. Read §0 first; it is the general
statement, and §§1–3 are the one instance of it that was actually run.

---

## 0. THE UNIFICATION — every bound this campaign can produce has the same shape

For any **centre** `D ⊆ {0..255}`, writing `k_D = Σ_{i∈D} 2^i`:

    k  =  k_D  +  sum_{i in S\D} 2^i  -  sum_{i in D\S} 2^i

so the residual target `T − fold(D)` is a **signed low-weight combination of ladder points**
exactly when `S` is close to `D`. Exhausting weight `M` against it proves

> ### `hamming_distance(S, D) > M`

* `D = ∅` → `w(k) > M`. **That is agent X's forward bound.**
* `D = {0..255}` → `256 − w(k) > M`. **That is this thread's complement bound.**

**The machinery is centre-agnostic: the cost is identical for every `D`.** So there is no
cleverness left to find in the search itself. The only thing that would make a particular `D` worth
testing is **a prior for `D`** — and the sole source of such a prior is how the instance was
constructed, which is closed by user instruction. That is where the remaining lever is, and it is
the user's to open.

### 0.1 Why the two centres already run are the only free ones — three independent proofs

**(i) XOR-affinity (this thread).** For a mask `m`, `k XOR m = k + m − 2·(k AND m)`. That is affine
in `k` — i.e. `(k XOR m)·G` is a fixed, computable transform of `k·G` — only when `k AND m` is
affine in `k`, which holds exactly for `m = 0` (gives `0`) and `m = 2^256 − 1` (gives `k`). **Of
2^256 masks, exactly two are usable.**

**(ii) Affine self-maps of `Z_N` (agent AB, independently).** AB reached the same uniqueness from
the group side without using the XOR identity.

**(iii) Signed-to-unsigned degeneracy (this thread, §4.2).** The general centre needs a *signed*
table. The two extremes are precisely the centres where the signed problem collapses to the
unsigned one — every `S` contains `∅`, and every `S` is contained in `{0..255}` — so they are the
only centres reachable with the table that already exists.

Three routes, one conclusion: **the upper bound has exactly one source.**

### 0.2 The earned result

> ## `10 ≤ w ≤ 246`
>
> Lower end: agent X, forward `w ≤ 9` exhausted, `2^33.1` candidates.
> Upper end: agent Y, complement `w' ≤ 9` exhausted, `2^33.1` candidates, 2,432 s.

Width 237 of 257 possible values. **Sharpening costs 42× per level and the cost is symmetric — there
is no cheaper side to push.** Against a `Binomial(256, ½)` null (mean 128, σ 8) the bracket excludes
about `2^-190` of the mass: **this is a certificate, not evidence.**

What it corrects is a claim in the record. `MINIMUM_COST_SEARCH.md` §7 said *"Upper bound: none
exists and none is obtainable from the instance."* **That was wrong.** One exists, it comes from
the same machinery as the lower bounds, and it scales at the same rate. The correction is worth more
than the number.

### 0.3 Running any centre you are given

`ycentre.py` is the driver. It emits engine inputs for a centre and needs **no target-specific
table** — the table depends only on the ladder. **Note the shared table `agentX_work/tbl4s.bin` +
`bm4.bin` was deleted by agent X at ~21:04 and must be rebuilt once before any further scan (§6.1);
after that, every centre reuses it.**

```
python3 ycentre.py 3,17,88,201     # explicit centre
python3 ycentre.py -               # empty centre   -> reproduces agent X's T   (CHECK: True)
python3 ycentre.py all             # full centre    -> reproduces T'            (CHECK: True)

./ymitm scan data_centre_<tag>_up.txt B ../agentX_work/tbl4s.bin ../agentX_work/bm4.bin rep.txt
   for B = 2,3,4,5   (covers distance 3..9; distance <= 4 is the direct probe in yedge.py)
```

**Read this caveat before quoting any centre result.** The engine here is **unsigned**, so it
certifies only the *one-sided* balls:

| file | target | what it certifies |
|---|---|---|
| `..._up.txt` | `T − fold(D)` | no `S ⊇ D` with `\|S \ D\| ≤ 9` |
| `..._down.txt` | `fold(D) − T` | no `S ⊆ D` with `\|D \ S\| ≤ 9` |

The **full two-sided Hamming ball** around a general `D` needs a signed table (`±L_i` on both
sides) and is **not** what this engine computes. For `D = ∅` and `D = {0..255}` the one-sided ball
*is* the full ball, which is why those two are unconditional.

**Validated end-to-end at a nontrivial centre** (`ycentre_test.py`): `D = {3,17,88,201}`,
planted `S = D ∪ {11,73,140,222,244}` (`|S| = 9`, distance 5); `T_fake − fold(D) == fold(E)` exactly,
and the unchanged engine recovered `E` with **all 10 splits FOUND**.

---

## 1. `T'` — the complement target, computed three ways

Curve, ladder and `T` re-derived from agent Q's raw instance-derived JSON (read-only) and
re-verified here; nothing inherited on trust (`ysetup.py` → `ydata.json`).

| check | result |
|---|---|
| `p == 2^256 − 2^32 − 977` (secp256k1 prime) | **True** |
| `N ==` secp256k1 group order | **True** |
| `a = 0`, 256/256 leaves on `y² = x³ + b` | **True** |
| `L_i == 2^i·G` by independent repeated doubling | **256/256, 0 bad** |
| `N·G == O`; `T` on curve; `N·T == O` | **True** |
| `G`, `T` identical to the values Q and X searched | **True** |

`A = (2^256 − 1)·G` by **three independent routes**: double-and-add with `(2^256−1) mod N`;
**folding all 256 ladder leaves**; double-and-add with the raw unreduced scalar. All three agree.

```
A.x  = 56062832069265557652411710355921883695174853966112376871108859166282132096323
T'.x = 34393883340176920870250405895813312293662019612542515124677151347541157631736
T'.y = 113211777249390963896039371927650543610689519859184504098294945786752299991906
```

| check | result |
|---|---|
| `T'` on curve; `N·T' == O`; `T + T' == A`; `T' ≠ T` | **True** |
| 12 random `S`: `fold(S)+fold(S̄)==A`, `k+k̄ == 2^256−1`, `fold==[k]G` | **12/12** |

**Independently reproduced by agent Z** from Z's own leaf extraction: `A` three ways agreeing and
matching mine, `T' = A − T` matching mine, `T + T' = A`, on curve, order `N`, complement identity
**20/20**. The construction does not rest on this thread alone.

## 2. The engine, and why no table was built

The MITM table holds the low 64 bits of `x(Σ_{i∈α} 2^i G)` for every `|α| ∈ [1..4]`.
**It depends only on the ladder, not on the target** — so agent X's `tbl4s.bin` (1.42 GB) and
`bm4.bin` (512 MB) worked unchanged. **No new large file was created.**
**Both files were deleted by agent X at ~21:04 and no longer exist — see §6.1 for the blast radius
(none: every result here predates the deletion) and for how to rebuild them.**
Revalidated independently before reuse (`ytable_check.py`):

| check | result |
|---|---|
| size == `8·(C(256,1)+…+C(256,4))` = 177,589,056 keys | **True** |
| sorted ascending (full scan of all 177.6 M keys) | **True** |
| bitmap exactly 2^32 bits | **True** |
| 60 random `\|α\|∈1..4` subsets, keys recomputed in Python bignum, present | **60/60** |
| 2,000 random 64-bit words present (negative control) | **0/2000** |

Engine: byte-identical copy of X's `xmitm.c` (`md5 1fec0554f59179c92e68eec7971193bd`) as `ymitm.c`,
recompiled. Field arithmetic re-cross-checked against Python bignum **on the complement base**
(`ycheck.py`): `x(L5+L7)`, `x(T'−L3)`, `y(T'−L3)` agree limb-for-limb.

## 3. RESULTS

### 3.1 Plant-test design — validate before any negative

`ycheckplant.py` requires **every** exact split of the planted set (within the scanned `i0` range)
to appear as a HIT — not merely "a hit occurred". A single hit can come from a carry-representation
and proves much less. Plant and real run differ **only** in the data file's first line.

| planted `S'` | `\|S'\|` | scan size (α size) | `i0` range | splits expected | found |
|---|---|---|---|---|---|
| `{11,73,140,201,244}` | **5** | 2 (α=3) | **full `[0,256)`** | 10 | **10 — PASS** |
| `{17,64,111,158,199,222,240,253}` | 8 | 4 (α=4) | `[17,18)` | 35 | **35 — PASS** |
| `{6,41,88,129,163,190,214,233,250}` | 9 | 5 (α=4) | `[6,7)` | 70 | **70 — PASS** |
| `{150,170,185,200,215,228,239,247,252}` | 9 | 5 (α=4) | `[150,151)` | 70 | **70 — PASS** |
| `{22,60,101,144,177,203,246}` | 7 | 5 (α=2) | `[22,23)` | 15 | **15 — PASS** |
| centre test, `E = {11,73,140,222,244}` | 5 | 2 (α=3) | full | 10 | **10 — PASS** |

The required weight-5 plant passed at full range; the expensive size-5 path passed at a low and a
high `i0` range and at two different α sizes.

### 3.2 Coverage and counts

Table `|α| ∈ [1..4]` × scan `|β| = b` covers `|S'| = b+1 … b+4`; `|S'| ≤ 4` closed by direct probe.

| `\|S'\|` | how | result |
|---|---|---|
| 0, 1 | `T'` is affine; full point compare vs all 256 `2^i·G` | **0/256** |
| 2–4 | `key(T')` probed in the table (truncation gives false positives only, so a miss is exact) | **absent** |
| 3–9 | scan sizes 2/3/4/5 | **0 hits** |

Verbatim `rep_comp.txt`:

```
DONE size=2 range=[0,256) n=32640      zero=0    0.0s     ( == C(256,2) )
DONE size=3 range=[0,256) n=2763520    zero=0    0.3s     ( == C(256,3) )
DONE size=4 range=[0,256) n=174792640  zero=0   23.5s     ( == C(256,4) )
DONE size=5 range=[0,256) n=8809549056 zero=0 2432.2s     ( == C(256,5) )
```

> **COMPLEMENT WEIGHT `w' ≤ 9` IS EXHAUSTED. 0 hits. 0 degenerate `dx=0` events at any size.**

### 3.3 Completion-evidence checklist — what a finished sweep looks like in writing

§3.2 was first written while size 5 was still running. That was premature. The claim is now earned,
and this is the evidence that earns it — **the pattern to reuse, not just the result**:

1. **status file terminal marker** — `yrun.status` ends `finished size 5 at 20:44:48` / `ALLDONE`;
2. **the process is gone, by PID not by name** — `yrun.pid` = 32218, `kill -0` → **DEAD**; no
   `ymitm scan data_comp` process exists;
3. **the engine's own total** — `DONE size=5 … n=8809549056`;
4. **independent re-derivation of that total** — the **252** `i0=… done` lines in `yrun_5.log` have
   candidate counts `C(255−i0,4)` summing to **exactly 8,809,549,056 = C(256,5)**;
5. **every gap explained** — the four `i0` with no line (252–255) are skipped by the engine's own
   `256−(i0+1) < SZ−1` guard and contribute `C(3,4)=C(2,4)=C(1,4)=C(0,4)=0`.

### 3.4 The citable conditional form, for any *unfinished* sweep

Derived by agent Z while auditing this thread mid-run. `i0` is the **smallest index of `β`**, and for
any `S'` with `|S'| ∈ [6,9]` one may always take `β` = the **five smallest elements of `S'`**
(leaving `α` of size 1–4, which the table covers). Therefore:

> completing `i0 ∈ [0, L)` proves exactly: **no complement set of size ≤ 9 contains an index `< L`**.

Uncovered fraction of weight-9 sets is `C(256−L,9)/C(256,9)`; recomputed here to check Z's figure:
`L = 96` → **1.3341%**, `L = 132` → **0.1257%**, `L = 256` → **0%**. **A conditional bound with a
stated fraction is citable; an unearned ✔ is not.**

## 4. Two consequences worth carrying

### 4.1 No signed-digit sweep on the ladder reaches the complement class

It is natural to assume signed digits subsume this family, since

    near-all-ones k = (2^256 - 1) - sum_j 2^{e_j} = 2^256 - 2^0 - sum_j 2^{e_j}

is only `w' + 2` signed terms — but **one exponent is 256**. Agent Z checked X's `xsigned.c`: its
digit loop runs `i < 256`, so the alphabet is `±2^e, e ≤ 255`, and under that alphabet **the minimum
signed weight of `2^256 − 1` is 42**. The class is outside any affordable depth. **The two searches
are complementary, not nested.** The fix is **agent AA's `±2^256` offsets**, which reach it without a
table rebuild; AA has been told.

### 4.2 The signed/unsigned degeneracy

See §0.1(iii) and the caveat table in §0.3. This is the operational reason a general centre is not
free even though its *cost* is identical: it needs a table this directory does not have.

## 5. Endomorphism orbit

`ylam.py` derives `beta` (cube root of 1 mod `p`) and `lambda` (cube root of 1 mod `N`) and verifies
`phi(x,y) = (beta·x, y) == [lambda]` on **8/8 random points**.

```
beta   = 55594575648329892869085402983802832744385952214688224221778511981742606582254
lambda = 37718080363155996902926221483475020450927657555482586988616620542887997980018
```

Twelve targets, all distinct, on the curve, of order `N`: `±T, ±phi(T), ±phi²(T)` (scalars
`±k, ±λk, ±λ²k`) and the complement `A − X` of each. A hit on any recovers `k` by a known
automorphism. `c_T == T'` and `T` is X's, so ten are new.

**`|S| ≤ 4` on all 12 (`yorbit_edge.py`): no hit; none equals any `2^i·G` (0/256 each).**

**`|S| ≤ 8` on the ten new targets** — `yorbit_run.sh 4`, scan sizes 2, 3, 4 =
177,588,800 candidates per target. See §6 for the per-target outcome.

**Coverage per unit cost, computed:**

| option | cost multiple | scalars covered | per unit cost |
|---|---|---|---|
| one more weight level, `w' ≤ 9 → ≤ 10` | `C(256,6)/C(256,5)` = **41.8×** | `C(256,10)/C(256,9)` = **24.7×** | **0.59** |
| the 12-target orbit at `w ≤ 9` | **12×** | **12×** | **1.00** |

Raw scalar count favours the orbit by 1.7×. **That is honest only about arithmetic:** under the one
hypothesis anyone holds — the designer picked a low-weight `k` — **ten of the twelve targets carry
~zero prior**, so the orbit is a hedge against a *different* designer choice (low-weight `λk`, or
low-weight `N−k`), not more evidence about the one at hand.

## 6. Files, and status at hand-off

| file | what |
|---|---|
| `ysetup.py` → `ydata.json` | re-derives + verifies curve, ladder, `T`, `A`, `T'` |
| **`ycentre.py`** | **the general driver — any centre `D`, no table rebuild** |
| `ycentre_test.py` | planted end-to-end test at a nontrivial centre |
| `ytable_check.py` | independent validation of X's reused table + bitmap |
| `ygen.py` → `data_comp.txt`, `data_fwd.txt` | engine inputs |
| `ymitm.c` / `ymitm` | byte-identical copy of X's engine, recompiled |
| `ycheck.py` | C-vs-Python arithmetic cross-check on the complement base |
| `yedge.py` | `\|S'\| = 0, 1` and the `\|S'\| ≤ 4` table probe |
| `yplant.py`, `ycheckplant.py` | planted-answer validation (all-splits criterion) |
| `yrun.sh`, `yrun.pid`, `rep_comp.txt`, `yrun_*.log` | the complement sweep |
| `ylam.py` → `yorbit.json`, `yorbit_edge.py`, `yorbit_run.sh` | the 12 orbit targets |

Nothing here is a `*.pkl`, so a container restart does not wipe it. **No large binaries were
created**; the 1.9 GB of tables lived in `agentX_work/` and **has since been deleted by agent X**.
To rebuild (needed for any further scan, ~3.4 GB peak):
`./ymitm table <any data file> 4 tbl4.bin`, `python3 ../agentX_work/xsort.py`,
`./ymitm bitmap <any data file> tbl4s.bin bm4.bin` — the table ignores the target.

| item | state |
|---|---|
| deliverable re-verified | 39,026/39,033, failing `[12231,12270,12350,14584,18673,22044,29125]` |
| `T'` construction | **verified** 3 ways; reproduced independently by Z |
| table reuse | **validated**; no new table built |
| plant tests | **6/6 PASS**, all-splits criterion |
| complement `w' ≤ 9` | **EXHAUSTED**, 0 hits, counts exact, evidence in §3.3 |
| **bracket** | **`10 ≤ w ≤ 246`** |
| orbit `\|S\| ≤ 4`, 12 targets | **done, no hit** |
| orbit `\|S\| ≤ 8` | **4 of 10 exhausted, no hit**; 6 never ran — shared table deleted mid-sweep (§6.1) |
| centre driver | **built and validated**; reproduces both known centres exactly |

### 6.1 Orbit `|S| ≤ 8` — per-target outcome, and a failure worth reading

**The shared table was deleted mid-sweep and my status file lied about it.** Both facts matter
more than the orbit result.

**What happened.** `agentX_work/tbl4s.bin` and `bm4.bin` were removed at ~21:03–21:04 (replaced by
`rt_*.bin`, mtime 21:04). **This was not agent X's error.** The coordinator has recorded that it
issued a "free what you no longer need first" instruction without checking who else read those
tables, after Z's audit had already reported them shared and AA's to be an identical multiset.
Recorded here so a successor reading the crash does not misattribute it. Every `ymitm scan` after that **segfaults instantly** — `mmap` on a
missing file returns `MAP_FAILED` and the code dereferences it unchecked. Confirmed by direct
reproduction: `data_comp.txt`, `data_negT.txt` and all six remaining targets now exit 139, with
**empty logs and no report file**.

**My `yorbit_run.sh` echoed `"$NM done"` after the inner loop without checking any exit code**, so
`yorbit.status` gained six lines claiming success for scans that never produced a byte — all
stamped the same second, `21:03:05`, which is what gave it away. That file is quarantined as
`yorbit.status.UNRELIABLE` with a `yorbit.status.README` beside it. **The script is fixed**: it now
checks the table exists before starting and tests every exit code. **Ground truth for this sweep is
`rep_orbit_<name>.txt` and nothing else** — the `DONE` lines carry candidate counts that are checked
against `C(256,b)` exactly by `yorbit_report.py`, which is why the false claims never reached the
table below.

**What this does NOT touch — the blast radius, established by timeline:**

| time | event |
|---|---|
| 19:58 | table validated (`ytable_check.py`: exact key count, sortedness, 60/60, 0/2000) |
| 19:59–20:03 | plant tests PASS — they *found their planted answers through this table* |
| **20:04:16 → 20:44:48** | **complement sweep, `DONE size=5 n=8809549056`** |
| 20:51:29 → 21:03:05 | orbit targets `negT, lamT, neglamT, lam2T`, all `n = 174792640 = C(256,4)` |
| ~21:03:05–21:04 | table deleted |
| 21:03:05 → | remaining six targets segfault |

The headline result completed **twenty minutes before** the deletion, and four orbit targets
completed after it started but before the file went away. A missing or truncated table cannot
produce a `DONE` line with an exact binomial count — it crashes, as we now see directly. Also,
`unlink` cannot truncate a live mapping: a run that successfully `open`ed the file holds a complete
view until it exits. **`10 ≤ w ≤ 246` is unaffected.**

**To finish the six**, the table would have to be rebuilt (`./ymitm table <any data file> 4
tbl4.bin`, `xsort.py`, `./ymitm bitmap …`, ~3.4 GB peak, ~15 min). **DECIDED — do not.** The
coordinator has ruled explicitly: the standing no-large-tables instruction holds, disk is at 69%,
and the orbit is a hedge carrying ~zero prior under the only hypothesis anyone holds (§5).

> **4 of 10 exhausted with 6 explicitly never-run is the deliberate final record**, and a better one
> than 10 of 10 bought by spending 4 GB on a hedge. **This question is settled; do not re-open it.**

*(table below generated by `yorbit_report.py` from the report files; a target counts as done only
with a `DONE` line for every one of sizes 2, 3, 4 and counts matching `C(256,b)` exactly)*

<!--ORBIT-->
| target | scalar | `|S| <= 8` | hits | degenerate | time |
|---|---|---|---|---|---|
| `negT` | `-k` | **exhausted** | 0 | 0 | 373.2 s |
| `lamT` | `lam*k` | **exhausted** | 0 | 0 | 530.9 s |
| `neglamT` | `-lam*k` | **exhausted** | 0 | 0 | 80.7 s |
| `lam2T` | `lam^2*k` | **exhausted** | 0 | 0 | 81.3 s |
| `neglam2T` | `-lam^2*k` | *incomplete* — sizes done: none | 0 | 0 | — |
| `c_negT` | `(2^256-1)+k` | *incomplete* — sizes done: none | 0 | 0 | — |
| `c_lamT` | `(2^256-1)-lam*k` | *incomplete* — sizes done: none | 0 | 0 | — |
| `c_neglamT` | `(2^256-1)+lam*k` | *incomplete* — sizes done: none | 0 | 0 | — |
| `c_lam2T` | `(2^256-1)-lam^2*k` | *incomplete* — sizes done: none | 0 | 0 | — |
| `c_neglam2T` | `(2^256-1)+lam^2*k` | *incomplete* — sizes done: none | 0 | 0 | — |

**4 of 10 targets complete.** Each complete target scanned 177588800 candidates (`C(256,2)+C(256,3)+C(256,4)`), counts checked exactly against the binomials; combined with the `|S| <= 4` probe of §5 this exhausts `|S| <= 8` on that target. Targets marked *incomplete* are **not** exhausted and must not be quoted as such.

### 6.2 Highest-value next experiment

**Not another level on either end.** 9 → 10 costs 42× on either side and moves `[10,246]` to
`[11,245]`, against a null centred at 128 — 2^58 short in both directions, symmetric cost.

**The experiment that changes something is the one that supplies a centre** (§0). The machinery is
built, validated, and centre-agnostic; the missing input is a prior for `D`, and its only source is
closed by user instruction.

Absent that, the best remaining spend is **agent AA's `±2^256`-offset signed-digit route** (§4.1):
the only cheap thing covering a class neither X's sweep nor mine reaches, and it needs no table
rebuild.
