# Agent AA — offset-shifted signed-digit MITM

**Angle:** every low-complexity search in this campaign tests one hypothesis — *`k` itself is
simple*. That is one point in a family. Since

> `k·G = T`  ⟺  `(k − c)·G = T − c·G`

a search for *"`k − c` is simple"* is the **identical machinery with one different point**. So a
list `C` of structured offsets multiplies the number of hypotheses tested by `|C|`.

**The structural fact that makes it cheap:** in a MITM the *table* side (all signed `a`-term
ladder sums) does **not** depend on the target. **The table is offset-independent; only the scan
side moves.** So build the table once, as big as affordable, and spend the per-offset budget on
the cheap half. That asymmetry is the whole experiment.

**Deliverable unchanged.** `solve_lab/best/new_instance_partial_39026.json` → **39,026 / 39,033**,
failing `[12231, 12270, 12350, 14584, 18673, 22044, 29125]`. I produced no better assignment.

**Scope, per the coordinator's correction (agent AB's Theorems A and B):** this work is a search
**for `k`**, not a route to an upper bound on `w`. A hit is a full solve; a miss is a weak,
citable exhaustion statement about that class and nothing more.

---

## 0. What was rebuilt and re-verified (nothing inherited on trust)

`aa_setup.py` re-derives and re-checks everything from `agentX_work/xdata.json` (read-only):

| check | result |
|---|---|
| `p == 2^256 − 2^32 − 977` | **True** |
| `N ==` secp256k1 group order | **True** |
| `L_i == 2^i·G` by independent repeated doubling | **256/256, 0 bad** |
| all 256 ladder points on the curve | **256/256** |
| `N·G == O`, `T` on curve, `N·T == O`, `lad[0] == G` | **all True** |

**Engine cross-check, the strongest one available:** my table generator, rewritten with sharded
output, produces a signed `a ≤ 3` table that is **bit-for-bit identical to agent X's independently
validated `stbls.bin`** — 11,119,616 keys, `(mine == X).all() → True`. So the shared half of the
machinery is not merely "believed correct", it reproduces another agent's artefact exactly.

---

## 1. The offset list `C` — 51 offsets, one line of rationale each

Ranked by **P(hit)** (per the coordinator), not by bound strength. `reach(c)` = minimum number of
signed terms with exponents ≤ 255 needed to write `±c mod N` — **the number of extra levels the
plain `c = 0` search would need to cover this offset for free.** Low reach ⇒ redundant.

### Tier 1 — a key-chooser plausibly typed this constant

| tag | `c` | reach | why a key might sit near it |
|---|---|---|---|
| `c0` | 0 | 0 | the plain case; baseline and the fleet's existing frontier |
| `a33` | `0x3333…33` | **128** | `(2^256−1)/5`, period-4 `0011` — a typed hex pattern |
| `a55` | `0x5555…55` | **111** | `(2^256−1)/3`, the densest period-2 constant |
| `aAA` | `0xAAAA…AA` | 103 | the other period-2 phase |
| `aCC` | `0xCCCC…CC` | 102 | the other period-4 phase |
| `inv10`,`inv5`,`inv3`,`inv7` | `1/d mod N` | 110/108/104/84 | `T = G/d` — a "nice" key in the *group*, not in binary |
| `cshift`,`bcurve` | instance constants | 92/90 | the generator may have reused a published constant as the scalar |
| `lam` | GLV λ | 83 | `k = λ + small` is `a + bλ` with `b = 1` and `a` of low weight — **not** covered by the campaign's `|a|,|b| < 2^21` exclusion |
| `lam2` | λ² | 84 | the third cube root of 1 |
| `drep` | `111…1` (77 digits) | 80 | decimal repunit — "nice" in base 10, structureless in base 2 |
| `d1e76`,`d1e38` | `10^76`, `10^38` | 68/33 | round decimal constants |
| `a11`,`a0F`,`a05`,`a03` | `0x1111…`,`0x0F0F…`,`0x0505…`,`0x0303…` | 64 | `(2^256−1)/15, /17, /51, /85` — period-4/8 patterns |
| `ones` | `2^256−1` | **42** | the complement: `k` with few ZERO bits. **agent Y owns this class** |
| `2p256`,`n2p256` | `±2^256` | 42 | `n2p256` **is `2N mod 2^256`**; these reach the exponent the 256-point ladder cannot hold (see §3) |
| `p2p256p1`,`n2p256m1` | `2^256+1`, `−(2^256−1)` | 41/42 | round-above-the-top constants |
| `halfN`,`inv2` | `(N∓1)/2` | 42 | half the group order / `1/2 mod N` |
| `pmodN` | `p mod N` | 42 | the field prime read as a scalar |
| `a01`,`aFF00` | `0x0101…`,`0x00FF00FF…` | 32 | `(2^256−1)/255`, `/257` |
| `a0001`,`aFFFF` | period 16/32 | 16 | `(2^256−1)/65535`, `/65537` |
| `a32_1`,`a32_F`,`a64_1` | period 32/64 | 8/8/4 | sparse repunits — **weakest members, near-redundant** |

### Tier 2 — negations (12): `n_a55, n_a33, n_a11, n_a0F, n_a01, n_aFF00, n_a0001, n_lam, n_inv3, n_d1e76, p2p257, n2p257`
`c` and `−c` are **different** offsets (`S(−c,m) = {N − c + small}` = "the group order minus a nice
constant"), so each needs its own scan. Ranked below the positives because "N minus a pattern" is a
less likely thing to type than the pattern.

### Tier 3 — controls, provably redundant (4): `p2_128, p2_255, p2_128p1, p2_255_1`
`2^128`, `2^255`, `2^128+1`, `2^255−1`. Included **only** to make the redundancy claim empirical as
well as proved. These are the coordinator's "`2^k ± 1` near 128" — see §2.

*Dropped as exact duplicates by the tooling, and recorded as results:* **`c = N` is bit-identical to
`c = 0`** (same base point), and `p2p256 == 2p256 mod N`.

---

## 2. The containment lattice — computed, not argued (`aa_lattice.py`)

Define the **signed-digit distance** `d(a,b) = w±(a−b)` = minimum number of terms `±2^e` with
distinct `e ≤ 255` summing to `a−b mod N` (minimised over both integer representatives in
`(−N,N)`). It is a metric, and

> **`S(c,m) = { k : d(k,c) ≤ m }` — every hypothesis in this family is a BALL, and the whole
> exercise is a packing problem.**

**(L0) Offsets only matter mod `N`.** `T − cG = T − (c+jN)G`. Two consequences the fleet should not
re-test: **`c = N` is exactly `c = 0`** (verified: identical base point), and **reduction is
invisible** — if the designer's integer `κ ∈ [N, 2^256)` has low signed weight, the plain `c = 0`
search already finds it, because the engine works on points and reduces for free. No offset is
needed for "the key was chosen ≥ N".

**(L1) `S(c,m) ⊆ S(0, m + reach(c))`.** (Minimal signed weight is subadditive: concatenate the two
representations and cancel/carry repeats, which never increases the count. Caveat: carries can push
an exponent from 255 to 256, which is exactly the off-by-one in §3.) Therefore:

- **`c = 2^e` for any `e ≤ 255`: reach 1 ⇒ `S(2^e,m) ⊆ S(0,m+1)`.** One extra level of the plain
  search subsumes **all 256 of them at once**. Testing `2^255` or `2^128` as an offset is
  **provably wasted budget.**
- **`c = 2^a ± 2^b` (so `2^128±1`, `2^255−1`, `2^k±1` for any `k ≤ 255`): reach 2 ⇒ `⊆ S(0,m+2)`.**
  Also provably wasted. *This is the coordinator's requested "`2^k ± 1` near 128" — it is dead by
  proof, and I ran two of them anyway as controls because they cost 30 s.*
- **`c = 2^256−1` has reach 42, not 2** — the ladder has no `2^256`, so the complement is genuinely
  non-redundant by 42 levels. **This is the same off-by-one agent AB flagged, detected independently
  here by the reach computation.**

**(L2) What `c = 0` already contains** (so nobody should test these as separate hypotheses):

| hypothesis | contained in | note |
|---|---|---|
| unsigned Hamming weight `w ≤ m` | `S(0,m)`, all `ε = +1` | signed at `m` **strictly dominates** unsigned at `w = m` |
| `k` with `≤ r` runs of ones | `S(0,2r)` | a run `b..a` is `2^{a+1} − 2^b`. **`m ≤ 7` ⇒ any `k` with ≤ 3 runs of ones**, e.g. `0x00FFFFF…F000…0FF`, unsigned weight ~40 |
| `k` close to a round constant `2^e` | `S(0,m)` | `2^e ± (m−1)` terms |
| short addition–subtraction chains | `S(0,m)` | length ≤ `m` from powers of two |
| any of the above on `κ ≥ N` before reduction | same | by (L0) |

**(L3) Offsets dominate basis extension.** A basis extended to exponent `E` covers, at depth `m`,
the `e = 256` branch with only `m−1` further terms; the offset `±2^256` at depth `m` covers it with
`m`. **And the offset costs nothing on the table side.** So I kept the 256-point basis and added
`±2^256, 2^256+1, −(2^256−1), ±2^257` as offsets — strictly stronger than rebuilding the table with
257 points, and no rebuild.

**(L4) `c` and `−c` are distinct** unless `2c` has weight ≤ `2m`. Hence tier 2.

**(L5) Measured packing — 1,275 pairs.** minimum distance 1; **1,222 pairs are at distance > 2m = 14
and are therefore PROVABLY DISJOINT at `m ≤ 7`**; 53 pairs overlap and are named in
`aa_lattice.py`'s output. Collapsing the overlap graph: **51 offsets ⇒ 34 well-separated clusters.**

Near-duplicates found (honest accounting of my own waste):

- **`lam2` and `n_lam` are distance 1** — because `λ² + λ + 1 ≡ 0 mod N`. The lattice computation
  rediscovered the GLV identity. One offset wasted.
- **`halfN` and `inv2` are distance 1** (`(N−1)/2` vs `(N+1)/2`). One offset wasted.
- **`ones` and `2p256` are distance 1**, `p2p256p1` distance 1 — so **agent Y's complement class and
  my `2^256` family are the same class at this depth.** Coordination note: Y should go *deeper* on
  `2^256−1` rather than sideways.
- The sparse repunits `a64_1, a32_1, a32_F, a0001, aFFFF` chain into `c0` at distances 4–16.

**(L6) The honest size of the whole thing.** `|S(c,7)| ≈ 2^50.6`; the union over 34 distinct classes
is `≈ 2^55.7`, i.e. **`2^−200.3` of the keyspace.** Against a uniformly random `k` this is worth
nothing, and no amount of offsets changes that. **The entire value is conditional on a designer
prior**, which is why the list is ranked by "what a person might type", not by cardinality.

---

## 3. Cost grid, and the frontier chosen

Candidate counts are exact (`C(256,b)·2^b`), timings measured on this 4-core box.

| side | `a` or `b` | candidates | cost |
|---|---|---|---|
| **table** (offset-independent, built ONCE) | `a ≤ 3` | 11,119,616 | 1.8 s, 89 MB |
| | **`a ≤ 4`** | **1,409,460,736** | **196 s build + 400 s sort + 24 s bitmap, 11.3 GB** |
| **scan** (per offset) | `b = 1` | 512 | <0.1 s |
| | `b = 2` | 130,560 | 0.4 s |
| | **`b = 3`** | **22,108,160** | **20–30 s** |
| | `b = 4` | 2,796,682,240 | ~45–60 min |

Reachable depth is `m = a + b`. The grid of what a `|C|`-offset sweep costs:

| table | scan | depth `m ≤` | per offset | 51 offsets |
|---|---|---|---|---|
| `a ≤ 3` (89 MB) | `b ≤ 3` | 6 | 4 s | 3.5 min |
| `a ≤ 3` | `b ≤ 4` | 7 | ~45 min | **38 h — unaffordable** |
| **`a ≤ 4` (11.3 GB)** | **`b ≤ 3`** | **7** | **~25 s** | **~25 min** |
| `a ≤ 4` | `b ≤ 4` | 8 | ~50 min | 42 h — pick a few |

**The chosen frontier: pay once for the deep table, then buy breadth.**
`a ≤ 4` moves `m ≤ 7` from *38 hours* to *25 minutes* for the whole sweep — a **~90× swing that
exists only because the expensive half is the offset-independent half.** That is the experiment.

**Breadth before depth, deliberately.** One extra level multiplies one class by ~62× in
*cardinality*; one extra offset multiplies the number of *structural families* by 1. Under a random
`k` neither matters (§L6). Under a designer prior the mass sits on **exact** structured constants
plus **very few** corrections — so a new family is worth more than a deeper ball on an old one.
Depth is then spent on the two or three highest-prior offsets.

---

## 4. Validation — planted answers, per offset class

**Requirement met before any negative was recorded.** `aa_setup.py` plants, for 23 offset classes, a
random 5-term signed `δ` and sets the fake target `T' = (c + δ)·G`; the pipeline then computes the
base `B = T' − c·G` **through the same code path as the real runs**, so offset bookkeeping is under
test, not bypassed. `aa_check_plant.py` predicts in independent Python the exact
`HIT <sz> <code> <s_last> <key>` line the engine must print, at **both** the `|α|=3/|β|=2` and
`|α|=2/|β|=3` splits, then decodes that line back to a scalar and re-verifies `k·G == T'` on the
curve.

**The sign-bookkeeping failure mode is real and this caught it:** my first predictor assumed the
table's "lowest-exponent sign forced to `+`" WLOG had to be undone in the scan indices. It does not —
only `x` is compared, so the `±` branch absorbs the global negation. **4 of 8 plants failed, exactly
the 4 whose lowest-exponent term was negative.** The engine was right and the predictor was wrong;
had I only checked "some hit appeared", the bug would have passed unnoticed in the predictor and I
would have had no independent model of the engine at all.

Planted `k` values have unsigned Hamming weights **40 to 188** — every one of them is invisible to
any plain-weight search ever run in this campaign, which is the coverage gain made concrete.

---

## 5. Files

| file | what |
|---|---|
| `aa_setup.py` → `aa_offsets.json`, `data/` | verification, offset list, per-offset base points, plants |
| `aa_signed.c` / `aa_signed` | the engine (X's `xsigned.c` + sharded table + bitmap mode) |
| `aa_field.h` | field arithmetic (copied verbatim from X's `xfield.h`) |
| `aa_sort.py` | per-shard sort |
| `aa_lattice.py` | the containment lattice / packing computation |
| `aa_check_plant.py` | planted-answer validation |
| `aa_run.sh`, `aa_wave1.sh` | sequential sweep drivers (ONE process at a time) |
| `tbl/t4.0 … t4.7` | **the a ≤ 4 signed table, 11.3 GB — coordinator: this is the large artefact** |
| `tbl/bm4.bin` | 512 MB prefilter bitmap (26.88 % of `2^32` bits set) |
| `runs/r_d_<tag>.txt`, `runs/r_p_<tag>.txt` | per-offset scan reports (real / plant) |

**Restart:** `gcc -O3 -march=native -fopenmp -o aa_signed aa_signed.c`; `python3 aa_setup.py`;
`./aa_signed table data/d_c0.txt 4 tbl/t4 4`; `python3 aa_sort.py tbl/t4`;
`./aa_signed bitmap tbl/t4 tbl/bm4.bin`; `./aa_wave1.sh`.
(numpy lives in `../agentX_work/pylib` — set `PYTHONPATH`.)

---

## 6. Results

### 6.1 Validation — **PASS for all 23 planted classes, at both splits**

`aa_check_plant.py` against the production `a ≤ 4` table:

> 23 offset classes × 2 splits = **46 predictions, 46 lines present, 46 decodes verified on the
> curve.** Planted `k` unsigned Hamming weights: **40 to 188.**

Classes validated: `c0, ones, a55, a01, inv3, lam, 2p256, d1e76, n2p256, n2p256m1, p2p256p1,
p2p257, n2p257, n_a55, n_lam, n_a0001, halfN, a0F, aFFFF, drep, a64_1, p2_128, p2_255_1` — i.e. at
least one member of **every structural family in `C`**, including both sign branches, both
exponent-256 branches, and both redundant-control families.

**The 8-shard lookup path was separately validated the same way**, because it is a different
lookup path and a negative from it would otherwise be unvalidated: plants `c0` and `lam` run through
all 8 shard passes each produced **5,595 hits — the identical count the monolithic-table run
produced** — and the predicted `HIT` line is present at both splits and decodes to the planted `k`.
Two lookup implementations, same answer to the last hit.

### 6.2 The sweep — counted, not claimed

**Every negative below is a negative of a search shown to find what it is looking for**, and every
"exhausted" cell below is a count the engine had to emit, checked against `C(256,b)·2^b`. The full
per-offset table is `FINAL_TABLE.txt`, produced by `aa_report.py`.

| | offsets | hits |
|---|---|---|
| **exhausted at `m ≤ 7`** | **41 of 51** | **0** |
| **exhausted at `m ≤ 6`** | **51 of 51** | **0** |
| never run | **0** | — |

**2,780,952,576 scan candidates counted. Degenerate `dx = 0` events: 0, at every offset and every
size. Expected false positives over the whole sweep: 0.212** — so a single `HIT` line would have
been overwhelmingly likely to be real. There were none.

The 10 offsets that stand at `m ≤ 6` rather than `m ≤ 7` are, deliberately, the least valuable ones:
`a32_F, a64_1` (reach 8 and 4 — the two most nearly redundant members), `n2p257`, `n_a01`,
`n_aFF00`, `n_a0001` (tier-2 negations of low-reach constants) and the four tier-3 controls that are
**provably** redundant anyway. Every high-reach offset, both `±2^256` branches, `c0`, `ones` and all
five high-reach negations are complete at `m ≤ 7`.

### 6.3 Two cost-model findings, both measured

**(a) A MITM table larger than the residual page cache is not a faster table, it is a disk-bound
one.** Thirty-two offsets ran at 15–30 s each while ~12 GB of cache was free. When the box's load
reached 26 and cache fell below the table's 11.3 GB, the running process showed
`utime=2113 stime=34946 majflt=690755`, climbing 6,000 major faults per 3 s — **94 % of its time in
the kernel servicing page faults**, turning a 20 s scan into a projected 50 min one. (PID obtained
by walking the process tree from the recorded wrapper PID; never by command-line matching.)

**(b) The fix is to shard the table and pass over it 8 times, and it is *faster* than the monolith
even when the monolith is warm.** Because `shard = key>>61` partitions the 64-bit key space, eight
passes against one 1.4 GB shard each have union exactly equal to one pass against the whole table.
Working set drops from 11.3 GB to **1.4 GB shard + 0.5 GB bitmap**, which stays resident under any
plausible pressure. Measured: **8.2 s per shard pass, 66 s per offset** — against 20–30 s warm and
~50 min thrashing for the monolith. It needed **no C change at all**: eight directories of symlinks,
one real shard and seven zero-length stubs, and the existing `st.st_size ? mmap : NULL` path already
does the right thing.

> **Recorded for the fleet: size your MITM table to the resident set, not to the disk. And when it
> does not fit, shard it — 8× the recursion is cheaper than 1× the page faults.**

### 6.4 Evidence discipline (coordinator's audit response)

`aa_shard.sh` originally had the defect agent AI found fleet-wide: unconditional `echo "SHARD$s"`,
engine stderr masked to `/dev/null`, and — worse than agent Y's case — **the resume guard keyed on
that same shell-written marker**, so a crashed shard would have been skipped forever. **Fixed
before the next batch**: exit code tested, stderr to `shardlogs/<tag>.b<b>.s<s>.err`, nothing but
the engine writes to the evidence file, and the resume guard keys on the engine's own
`DONE … n=<count>` line checked against `C(256,b)·2^b`.

**Item 3, the check rather than the expectation:** re-verified the already-marked shards against
engine output alone — **234 shell markers, 234 engine `DONE` lines carrying an exact closed-form
count; markers never outran evidence; 9 of 10 offsets cleared 24/24, the tenth (`n2p257`) was
mid-run at 18/24 and is reported as `m ≤ 6`.** No damage, as AI predicted, but the count is here.

`aa_report.py` **ignores every shell marker by construction** and prints `NEVER RUN` rather than
omitting a row. The build was also re-done without the `gcc … | head && echo BUILD_OK` mask
(`gcc exit=0`, tested directly).

### 6.5 Negatives, one line per offset — and each is weak on its own

For every offset `c` in the list, the statement earned is exactly:

> **`k − c` has signed-digit weight `> m`** (with `m = 7` for the 32 offsets in the first sweep,
> `m = 6` for the remaining 19), where the digits are `±2^e` with `0 ≤ e ≤ 255`.

Individually each of these is nearly worthless: `|S(c,7)| ≈ 2^50.6` against a `2^256` keyspace, and
`Pr[hit] ≈ 2^−205` under a uniform `k`. **The value is the union of classes covered, and even that
union is `2^−200.3` of the keyspace.** What the sweep actually buys is the elimination of an entire
*shape* of designer hypothesis — "the key is a well-known constant, give or take a handful of bits" —
across 34 structurally distinct constants at once, rather than one more level on the single
hypothesis the fleet has been testing all along.

Specifically now excluded (all with `≥ 6` signed terms of slack):

- `k` is **not** any repunit-type constant `(2^256−1)/{3,5,15,17,51,85,255,257,65535,65537,2^32±1,2^64−1}`
  or its complement/negation, plus up to 6–7 signed corrections;
- `k` is **not** `1/2, 1/3, 1/5, 1/7, 1/10 mod N`, `(N±1)/2`, `λ`, `λ²`, `N−λ`, `p mod N`, the
  instance's `b` or `c_shift` read as scalars, `10^38`, `10^76`, or the 77-digit decimal repunit,
  plus corrections;
- `k` is **not** `±2^256`, `±(2^256−1)`, `2^256+1`, `±2^257` plus corrections — **the region the
  256-point ladder structurally cannot represent at all**, and the reason that region needed offsets
  rather than a deeper plain search;
- and by (L1), **not** `2^e ± 2^f` plus 5 corrections for any `e, f ≤ 255`, which the two
  redundant-control offsets confirmed empirically as well as by proof.


---

## 7. Check-in 117 follow-up — the approved `m ≤ 8` run on `c0`

### 7.1 Duplicate pruning (coordinator's ruling — agent Y's thread is closed, so I took it)

My own lattice showed the fleet paying more than once for one class. Pruned unilaterally:

| cluster | pairwise distance | **kept** | dropped |
|---|---|---|---|
| `ones = 2^256−1`, `2p256 = 2^256`, `p2p256p1 = 2^256+1` | 1 | **`2p256`** | `ones`, `p2p256p1` |
| `n2p256 = −2^256`, `n2p256m1 = −(2^256−1)` | 1 | **`n2p256`** | `n2p256m1` |
| `halfN = (N−1)/2`, `inv2 = (N+1)/2` | 1 | **`inv2`** | `halfN` |
| `lam2 = λ²`, `n_lam = N−λ` (because `λ²+λ+1 ≡ 0`) | 1 | **`lam`** family: keep `lam2` | `n_lam` |

`2p256` and `n2p256` are kept because they are the *pure* exponent-256 terms — the exact centres
of their clusters, with the others sitting at distance 1. **The precise pruning statement:
running the representative at depth `m` subsumes every dropped member at depth `m−1`.** All of
them are already complete at `m ≤ 7`, so nothing is lost; and it means any future depth on
`2p256` automatically re-covers `ones` one level behind, which is where agent Y's class now lives.

**Per the ruling, no further `2^e` / `2^a ± 2^b` control offsets will be run.** Two were run, they
behaved exactly as `reach = 1, 2` predicts, and that is the end of it — the rest is a theorem.

### 7.2 Validation FIRST, at the awkward cases (`aa_plant8.py`)

`m = 8` is reachable by **exactly one split**, `|α| = 4 / |β| = 4`, so a plant must be found at
`b = 4` or not at all. Four cases chosen to be the ones that break if bookkeeping is subtly wrong:

| case | what it stresses | unsigned wt of `k` | result |
|---|---|---|---|
| `allneg` | **every sign negative** — table's lowest-exponent term negative, so the match must come out on the `−` branch. *This is the case that caught my own predictor.* | 191 | **predicted line found, bit for bit** |
| `edges` | exponents **0 and 255** both present — both ends of the ladder | 142 | **found** |
| `adjacent` | adjacent pairs `(e, e+1)` in **both** halves — stresses `nxt(s) = ((s>>1)+1)<<1`, which is what forbids repeated exponents | 119 | **found** (1,175 hits: the expected carry-representation multiplicity) |
| `tight` | all eight exponents in a narrow **high** window — stresses the recursion cut `nxt(s) > 512−2·(SZ−1−depth)`, the one place an over-eager prune silently loses answers | 10 | **found** |

All four verified across **8/8 shard passes**, predicted `(sz, code, s_last, key)` computed in
independent Python before the engine ran, and `k·G == T'` re-checked on the curve.
`allneg` came out on branch `−1` exactly as predicted.

### 7.3 The run

`aa_deep8.sh`: outer loop = **16 balanced `s0` chunks**, inner loop = **8 table shards**. All eight
shards of a chunk finish before the next chunk starts, so *"chunks 0..j complete"* is an **exact
fraction of the m = 8 candidate space**, not a vague partial. Total 2,796,682,240 candidates × 8
shard passes = **22.4 × 10⁹ candidate-evaluations**. Evidence-keyed resume, exit codes tested,
stderr to `shardlogs8/`.

**Status at hand-off — IN FLIGHT, and reported as an exact fraction, never as a claim.** First
unit closed clean: `DONE data/d_c0.txt sz=4 range=[0,9) n=192557240 zero=0 320.7s` — the count is
exactly the closed form for that range. **1 of 128 shard-units complete, 0 hits, 0 degenerate
`dx = 0` events.** Measured 600 k candidate-evals/s under load ~20 (against 2.7 M/s uncontended),
so the projection is **~43 min per chunk, ~10.4 h for the full `m ≤ 8`**. Coverage is reported by
`aa_prog8.py`, which counts a chunk **only** when all eight shard passes carry the engine's own
exhaustive `DONE` line, and lists partially-run chunks separately without counting them.
Engine PID recorded in `deep8_pid.txt`, located by **open-file ownership** (`/proc/*/fd`), not by
command-line matching.

*(A note against my own earlier confusion: the `.pid` files that "vanished" were not swept by the
environment — `cd X && setsid … &` backgrounds the `cd` too, so `echo $! > f` wrote the file in the
original working directory. My bug, found and stated.)*

### 7.4 What it buys, stated against a model rather than as a bound on `w`

Per agent AC's standard, and **not** claimed as progress on bounding `w`, which it is not:

- **Under uniformity:** `|S(0,8)| = 2^56.56` (a 62.4× ball, from `2^50.60`), so
  `P(hit) = 2^−199.4`. That is the whole story against a random key and it is nothing.
- **Against the stated alternative `H_W`** — *designer picks the signed-digit weight uniform on
  `{1..W}`* — a miss at `m ≤ M` removes `M/W` of `H_W`, and the likelihood ratio against `H_W`
  is `1/(1−M/W)`:

| `W` | `m ≤ 7` (already had) | **`m ≤ 8` (this run)** | marginal gain |
|---|---|---|---|
| 8 | 8.00× | **killed outright** | kills it |
| 10 | 3.33× | **5.00×** | 1.50× |
| 12 | 2.40× | **3.00×** | 1.25× |
| 16 | 1.78× | **2.00×** | 1.12× |
| 20 | 1.54× | **1.67×** | 1.08× |
| 30 | 1.30× | 1.36× | 1.05× |
| 60 | 1.13× | 1.15× | 1.02× |

> **So the honest value of this run is: it kills `H_8` outright, and beyond `W ≈ 20` it moves
> nothing worth reporting.** The `2^50.6 → 2^56.6` sixty-fold growth in covered keys converts to a
> 1.08× nudge at `W = 20`. That is the exchange rate, and it is why depth stops being worth buying
> not far past here.

*(Standing prohibition respected: the exchange rate is computed against a stated model of the
designer's choice. Nothing here investigates how the instance was actually produced.)*
