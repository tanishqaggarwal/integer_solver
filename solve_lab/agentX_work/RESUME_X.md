# Agent X — low-Hamming-weight meet-in-the-middle on the fold scalar

**Angle:** self-contained search with a clean pass/fail. A hit ends the whole problem
(39,033/39,033 via T's integer lift); a miss is a citable exhaustion bound on `k`.

**Deliverable unchanged and re-verified at start:**
`solve_lab/best/new_instance_partial_39026.json` → **39,026/39,033**, failing
`[12231, 12270, 12350, 14584, 18673, 22044, 29125]`. I produced no better assignment.

---

## 1. WHAT WAS REBUILT AND RE-VERIFIED (nothing inherited on trust)

`xsetup.py` re-derives everything from agent Q's raw instance-derived files (read-only) and
re-checks it:

| check | result |
|---|---|
| `p == 2^256 − 2^32 − 977` | **True** — this is the secp256k1 prime |
| `N == 0xFFFF…D0364141` | **True** — this is the secp256k1 group order |
| 256 leaf points on the cubic `y² = x³ + b`, `a = 0` | **256/256** |
| `L_i == 2^i·G` by independent repeated doubling from `G = leaf(x2779)` | **256/256, 0 bad** |
| `N·G == O` | **True** |
| `T = (C1 + K/3, C2) mod p` on the curve, `N·T == O` | **True** |
| `G`, `T` equal the values Q's six search programs used | **True** |

So the curve is **isomorphic to secp256k1** (same `p`, same `N`; `b = u⁶·7`). That is an
independent confirmation of "no exploitable structure": prime order, non-anomalous, no small
embedding degree, endomorphism worth √3.

Artifacts: `xdata.json` (curve, G, T, all 256 ladder points), `data_real.txt` (engine input).

## 2. THE ENGINE

Q's Python MITM ran at ~24k candidates/s (`wt7.log`: 59.9M in ~2,500 s). At that rate `w ≤ 9`
is ~100 core-hours. I rewrote the inner loop in C: `xmitm.c`, 256-bit field arithmetic
specialised to `p = 2^256 − c` (pseudo-Mersenne fold, `__int128` limbs), **affine** addition with
**batched inversion** (5 field mults per candidate; only `x3` is computed at the last level),
Gray-style nested recursion so each candidate costs one addition, OpenMP over the outermost index.

Measured: **~250 ns/candidate** single-threaded for table build; **~1.23 µs/candidate wall**
for the scan under contention with the rest of the fleet on this 4-core box. That is **~50×**
Q's rate.

Lookup: table = **sorted 64-bit truncated x-coordinates only** (no codes — a hit is re-decoded
by rescanning, since hits are ~0), 1.42 GB, plus a **2^32-bit (512 MB) prefilter bitmap**
(`bm4.bin`) so 96% of probes cost one memory access. Expected false positives over the entire
`w ≤ 9` sweep: **0.085**.

## 3. VALIDATION — the machinery passes the planted-answer test

**Field/point arithmetic vs Python bignum:** `x(L5+L7)`, `x(T−L3)`, `y(T−L3)` agree limb-for-limb.

**Table generation vs Python group law:** size-1 block **256/256** exact, size-2 block
**32,640/32,640** exact, size-3 block **400 random samples, 0 mismatches**.

**Planted targets** (`xplant.py` sets `T' = Σ_{i∈S} 2^i·G` for a chosen `S`, scan range restricted
to the planted `i0`):

| planted `S` | \|S\| | scan size | found? |
|---|---|---|---|
| `{3,41,97,180, 200,214,229,241,250}` | 9 | 5 | **yes, uniquely** |
| `{7,50,133, 200,214,229,241}` | 7 | 4 | **yes** |
| `{7,50,133,199, 200,214,229,241}` | 8 | 4 | **yes** |
| `{11, 200,214,229,241,250}` | 6 | 5 | **yes** |
| `{200,214}` | 2 | 2 | **yes** |

`xdecode.py` decoded the weight-9 hit back to **α = [3,41,97,180], β = [200,214,229,241,250]** —
exactly the planted set. **The negatives below are therefore negatives of a search that has been
shown to find what it is looking for.**

*(The extra hits on the \|S\| ≤ 6 plants are not false positives: they are genuine
carry-representations, e.g. `2^250 = 2·2^249`, giving a different (α,β) with the same `k`.
`xdecode.py` reconstructs `k = Σα + Σβ` and re-verifies `k·G == T` exactly, so these are
handled, not miscounted.)*

## 4. RESULTS

Coverage design: table = every subset of size **1..4** (177,589,056 keys); scan = every subset
of size **1..5**. Every `S` with `2 ≤ |S| ≤ 9` splits as `α ⊎ β` with `1 ≤ |α| ≤ 4`,
`1 ≤ |β| ≤ 5`, so the pair (table, scan) covers it. `|S| = 0` and `|S| = 1` checked directly
(`xedge.py`: `T` is not `O` and is not any `2^i·G`). Degenerate `dx = 0` events are trapped and
reported, not silently skipped — **0 occurred**.

| weight | candidates | status | time |
|---|---|---|---|
| \|S\| = 0, 1 | 257 | **exhausted, no solution** | — |
| \|S\| ≤ 6 | — | **exhausted, no solution** (re-done, subsumed) | — |
| \|S\| ≤ 7 | — | **exhausted, no solution** (Q's 33.7% partial now CLOSED) | — |
| \|S\| ≤ 8 | scan sizes 2+3+4 = 177,588,800 | **exhausted, no solution** | 219 s |
| \|S\| ≤ 9 | scan size 5 = 8,809,549,056 | **RUNNING** | est. ~3 h |

Zero-events 0, hits 0 at every completed size.

## 5. FILES

| file | what |
|---|---|
| `xsetup.py` → `xdata.json` | re-derives + verifies curve, G, T, ladder |
| `xgen.py` → `data_real.txt` | engine input |
| `xmitm.c` / `xmitm` | the C engine (`table` / `scan` / `bitmap` / `selftest`) |
| `xcheck1.py`, `xcheck2.py` | arithmetic + table cross-checks against Python |
| `xplant.py`, `xdecode.py` | planted-answer test and hit decoding |
| `xedge.py` | \|S\| = 0, 1 and the `b = 0` edge case |
| `xsort.py` | sorts `tbl4.bin` → `tbl4s.bin` |
| `tbl4.bin` (1.42 GB, build order) / `tbl4s.bin` (sorted) / `bm4.bin` (512 MB bitmap) | the table |
| `runAB.sh`, `scanAB.pid`, `scanAB.log`, `rep_real.txt` | the running sweep |
| `K_CONSTRAINTS.md` | the verified constraint catalogue on `k` |

**To resume after a restart:** `gcc -O3 -march=native -fopenmp -o xmitm xmitm.c`;
`python3 xsetup.py`; `python3 xgen.py`; `./xmitm table data_real.txt 4 tbl4.bin`;
`python3 xsort.py`; `./xmitm bitmap data_real.txt tbl4s.bin bm4.bin`; then
`./xmitm scan data_real.txt 5 tbl4s.bin bm4.bin rep_real.txt [i0lo] [i0hi]`.
The `[i0lo,i0hi)` range makes the size-5 sweep restartable in pieces — `rep_real.txt` records a
`DONE size=5 range=[lo,hi)` line per completed piece.
(`pylib/` holds locally-installed numpy + gmpy2; the restart had wiped both.)

---

## 6. WORK ADDED ALONGSIDE THE SWEEP (coordinator's second task)

### 6.1 `K_CONSTRAINTS.md` — the verified constraint catalogue on `k`
Every row re-derived or re-executed here; rows only code-audited are marked as such. Includes the
standing ruling (Q's searches **do** have instance-level standing for the negative direction, with
the one residual risk named: F's parse, not the modulus) and the per-bit analysis.

### 6.2 NEW RESULT — Q's slot-collision caveat is vacuous for every `|S| ≤ 42`
Q left this open in its §4 ("needs the particular scalar"). It does not. Two children of a stage
slot coincide iff `Σ_{S1} 2^i − Σ_{S2} 2^i = ±N` with `S1, S2` disjoint subsets of the ON-set — that
is a signed-binary representation of `N` with `|S1|+|S2|` nonzero digits, and the minimum over all
such representations is the **NAF weight of `N`**, which I computed to be **43**
(`weight(N) = 192` in plain binary; the NAF reconstructs `N` exactly). So for any `|S| ≤ 42` no slot
ever sees two equal live children and the degenerate branch — where a gadget imposes **nothing**
because `dx = dy = 0` — never opens. This makes the forward implication
*(satisfying assignment ⟹ `k_S·G = T`)* airtight across the entire low-weight regime that this
fleet works in, including T's integer lifts at `|S| = 2,3,5,6,7,8,17`.

### 6.3 STRENGTHENED — BSGS bound pushed from 2⁴⁴ to 2⁵²
`xbsgs.c`, same field code as `xmitm.c` (`xfield.h`), 2²⁶ baby steps × 2²⁶ giant steps in W = 512
parallel chains with batched inversion. `smul` cross-checked against Python at 5 scalars, and a
**planted `k₀ = 5·2²⁶ + 1234567` was recovered at exactly `i = 5`**.
Result: `+T` → **0 candidates**. `−T` → see `bsgs.log`.
So `k > 2⁵²` and `N − k > 2⁵²`. Baby table: `babys.bin` (537 MB, sorted).

### 6.4 Per-bit measurement on the instance (`xperbit.py`)
Over the 256 leaf-selector wires in `EQUATIONS.txt`: **0 / 256** appear in an equation mentioning no
other wire, so **no selector is pinned**. Their footprints do vary (77–185 occurrences across 30–51
equations) — that is *where a leaf sits in the fold tree*, identical for every candidate `S`, and
carries no information about whether the bit is ON.

### 6.5 Note on scheduling
The box is heavily oversubscribed by the rest of the fleet (load 22–28 on 4 cores); the size-5 scan
gets ~0.4 cores despite 8 OMP threads. It is **not** I/O bound (`read_bytes = 0`, RSS 1.9 GB fully
resident) — purely CPU-starved. `tbl4.bin` (unsorted, 1.4 GB) and `baby.bin` were deleted to relieve
page-cache pressure; both are regenerable in under a minute.
