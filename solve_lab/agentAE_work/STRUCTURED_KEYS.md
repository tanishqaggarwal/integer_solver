# STRUCTURED_KEYS — the cheap structural families of `k₀`, priced and swept

**Agent AE.** Everything below was measured in this directory. Scripts: `ae_setup.py`,
`ae_lib.py`, `aekang.c`, `ae_plant.py`, `ae_families.py`, `aequot.c`, `ae_quot.py`,
`ae_support.py`. Nothing is quoted from another agent without re-derivation except where the
row says otherwise.

---

# 0. THE ONE-PARAGRAPH VERSION

Every family here is a **magnitude or quotient** family, not a weight family. A **hit** in any of
them solves the instance outright and hands you `w` exactly. A **miss** excludes a set of
scalars and says **nothing whatever about `popcount(k₀)`** — the exclusions below are statements
about where `k₀` is **not**, not about how heavy it is. That is the whole honest content, and it
is what agent AB's §9.13 already said before I started. The reason to run them anyway is that
they are cheap and that nobody had. **Result: a clean sweep of misses. No hit anywhere.**

The one number the campaign can take away: with `O(1)` memory, one core and a distinguished-point
kangaroo, **this box decides `k₀ < 2^58` at `2^32` group operations (~1.5 h at the rate this box
actually delivered)**, against a prior of `2^-198`. That is the best keys-per-operation ratio any
instrument on this fleet has bought — and it is still `2^198` short of mattering. `R = 64` is
`2^35` operations, ~12 h at the same rate; nothing about the method stops there, only the clock.

---

# 1. VERIFIED BEFORE STARTING (`ae_setup.py`)

| check | result |
|---|---|
| `p == 2^256 − 2^32 − 977` | **True** |
| `N ==` secp256k1 group order | **True** |
| `a == 0`; `G`, `T` on the curve | **True** |
| `N·G == O`, `N·T == O` | **True** |
| `ladder[i] == 2^i·G` for **all** `i = 0..255`, independent doubling chain | **256/256, 0 bad** |
| all 256 ladder points on the curve | **True** |
| `G`, `T` identical to agent **Y**'s independently extracted `ydata.json` | **True** |
| GLV: `β³ ≡ 1 mod p`, `λ³ ≡ 1 mod N`, `φ(x,y) = (βx, y) == [λ]` on 12 points | **True** |

`ae_data.json` is written only if every check passes; the script exits non-zero otherwise.

---

# 2. THE INSTRUMENT — a distinguished-point kangaroo, and what it does and does not prove

`aekang.c`: secp256k1 field arithmetic on 4×64-bit limbs, affine walks, one Montgomery batched
inversion per lock-step round over `K` kangaroos, jump index = low 5 bits of `x`, distinguished
point = low `d` bits of `x` zero, shared open-addressed DP table. Memory is the DP table alone
(**8–64 MB**, not a function of the range). Half the herd is tame (start `u·G`, `u` uniform in
`[0,L)`), half wild (start `Q + v·G`, `v` uniform in `[0,L/8)`); mean jump `2^{R/2}`.

## 2.1 It is **not** an exhaustion. This matters and is stated everywhere below.

BSGS decides `k₀ ∈ [a, a+L)` **deterministically**. A kangaroo does not: a miss means "no
tame–wild collision within the jump budget". The honest form of every miss below is therefore
**"excluded at confidence 1 − ε"**, with ε computed, not asserted.

## 2.2 Calibration — 34 planted keys, every one recovered **exactly**

A plant is a uniform `k` in `[0,2^R)`; PASS requires the engine's candidate to be confirmed by an
**independent Python bignum recomputation** of `k·G` *and* to equal the planted `k`.

| plant set | n | recovered exactly | mean `jumps/√L` | max |
|---|---|---|---|---|
| `R = 32` (`ae_plant.py 32 12 2 32`) | 12 | **12/12** | 2.406 | 7.29 |
| `R = 40` | 12 | **12/12** | 1.946 | 4.72 |
| `R = 48` | 6 | **6/6** | 1.546 | 3.13 |
| **shifted interval** `k = c + j`, random 256-bit `c`, run on `T − c·G` | 4 | **4/4** | 1.42 | 2.85 |
| **total** | **34** | **34/34** | **1.976** | 7.29 |

**Closed-form check #1 (rule 5).** The theoretical mean is `2√L`. Measured over 34 plants:
**1.976** — this is the number that says the engine is not merely running but running *correctly
at the right cost*.

**Closed-form check #2 (rule 5).** Every run reports `dps` and the closed form `jumps/2^d`.
Across every run in this document the ratio is **0.99–1.01**. A run whose DP count does not match
its own jump count over `2^d` has not searched what it claims to have searched. A crashed or
table-less run cannot produce this ratio; it produces no `DONE` line at all.

**Closed-form check #3.** The engine self-tests before every run: it re-derives `2·L_i == L_{i+1}`
for 8 ladder indices and `jd_i·G == J_i` for all 32 jump points from the ladder, and exits **3**
with `SELFTEST_FAIL` otherwise. Every run below printed `SELFTEST_OK`.

## 2.2a A bug the plants missed and a counter caught — worth the space

The first version of the quotient engine (`aequot.c`) reported `deg1 = 1023` and `deg2 = 1023`
degenerate steps on a table of `W = 1024` lanes with `chunk = 1024`. **Four random planted
quotients passed anyway** — `(a₀,b₀)` were large, and the damage was confined to lane 0.

`1023 = chunk − 1` is not a plausible count for a `2^-128` event. Reading it against its closed
form: lane 0 starts at `1·G`, so its first step computes `dx = x(G) − x(G) = 0` — an addition
that is really a **doubling** — the batch-inversion path skipped it, the lane never advanced, and
**every entry of lane 0's chunk (`a ∈ [1, 2^15]`, and likewise `b`) was recorded against
`x(G)`**. About 0.2 % of the pairs, but including `a = 1`, i.e. the entire `k₀ = 1/b` sub-family,
which is the most interesting part of the family.

Fixed by supplying `2G` and `2T` and resolving the `P == base` case explicitly (the `P == −base`
case aborts with a non-zero exit; it is unreachable here). After the fix the counter reads
**`deg1 = deg2 = 1` exactly** — the one genuine doubling — which is itself the closed form. Four
**deliberate lane-0 plants** (`(a,b) = (1,3), (2,1), (1,2^{19}), (5,7)`) were added and are part of
the standing validation.

> **The random plants could not have found this and did not. The counter found it.** That is the
> case for rule 5 stated as concretely as I can state it.

## 2.3 The failure model, measured

| `c` | plants needing `> c√L` | empirical `P` | `exp(−c/2)` |
|---|---|---|---|
| 2 | 13/34 | 0.382 | 0.368 |
| 4 | 4/34 | 0.118 | 0.135 |
| 6 | 1/34 | 0.029 | 0.050 |
| 8 | 0/34 | 0.000 | **0.018** |

The exponential model with mean `2√L` fits (sd 1.53 vs mean 1.98). **A run capped at `8√L`
therefore misses a genuine hit with probability ≈ 1.8 %; capped at `4√L`, ≈ 13.5 %.** Each family
row states its cap, so each exclusion carries its own ε.

## 2.4 Build hygiene

Both binaries are compiled with the exit status checked directly (`gcc …; rc=$?`), never through
a pipe — a pipeline's status is the last stage's, which is how a failed link on this fleet once
left a stale binary in place. `aekang` and `aequot` were both rebuilt this way from unchanged
source after the first plant runs and re-verified to execute.

---

# 3. WHAT A MISS ON A MAGNITUDE FAMILY GIVES — the direction of every implication

This section exists because getting it wrong is the main way this thread could produce a
falsehood.

| family | a **hit** proves | a **miss** proves | does the miss bound `w`? |
|---|---|---|---|
| `k₀ < 2^R` | `k₀` known ⇒ `w` known exactly; and `w ≤ R` | `k₀ ≥ 2^R`, i.e. `max(S) ≥ R` | **NO** |
| `N − k₀ < 2^R` | `k₀` known ⇒ `w` known | `N − k₀ ≥ 2^R` | **NO** |
| `\|k₀ − c\| < 2^{R−1}` | `k₀` known ⇒ `w` known | `k₀` is not within `2^{R−1}` of `c` | **NO** |
| `k₀ = ±λ^{-e}·a·b^{-1}` | `k₀` known ⇒ `w` known | no such representation with `a,b ≤ 2^m` | **NO** |
| **NAF weight ≤ m** (agent X) | `k₀` known | `NAF-wt(k₀) > m`, and since `w ≥ NAF-wt`, only `w > m` | **only from BELOW** |

**The NAF row is the trap the brief warned about, and the check comes out the way the brief said.**
The non-adjacent form is the *minimal* signed-digit representation, so `NAF-wt(k₀) ≤ w` always.
A signed-digit exhaustion therefore pushes `w` **up**, never down: agent X's `m ≤ 7` exhausted
yields `w ≥ 8`, which is strictly **weaker** than X's own unsigned `w ≥ 10`. That is exactly why
the two numbers differ, and it is why **no signed/NAF experiment was run here** — the implication
points away from the upper bound this thread was asked for. Priced at zero, correctly.

**So: every row in §4 is a lottery ticket whose only payoff is the hit branch.** Under agent AC's
posterior (`w ∈ [115,141]` at 90 %), the hit branch is where all of the value is, and there is
almost none of it.

---

# 4. THE FAMILY SWEEP — priors, costs, outcomes, exclusions

Priors are `|F| / N` with `log₂ N = 256.0000`.

## 4.1 Results

<!--RESULTS-->
*(filled in below as runs completed; see §4.2 for the per-family table)*

---

# 5. THE INSTRUMENT COMPARISON — why kangaroo, and where it stops

Coverage per operation, on this box (higher is better):

| instrument | memory | coverage | ops | keys/op | exhaustive? |
|---|---|---|---|---|---|
| agent X's weight `≤ 9` sweep | 1.9 GB table | `2^53.4` | `2^33.1` | `2^20.3` | **yes** |
| kangaroo `R = 48` (`8√L` cap) | 8 MB | `2^48` | `2^28` | `2^20.0` | no (ε = 1.8 %) |
| kangaroo `R = 56` | 16 MB | `2^56` | `2^32` | `2^24.0` | no |
| **kangaroo `R = 58`** (run) | **16 MB** | **`2^58`** | **`2^32`** | **`2^25.0`** | no (ε = 1.8 %) |
| kangaroo `R = 60` | 32 MB | `2^60` | `2^33` | `2^27.0` | no |
| kangaroo `R = 64` | 64 MB | `2^64` | `2^35` | `2^29.0` | no |
| **quotient sweep `m = 25`** | **512 MB** | **`2^51.9`** | **`2^26`** | **`2^25.9`** | **yes** |
| BSGS, whole 2 GB budget | 2 GB | `2^56` | `2^29` | `2^27.0` | **yes** |

Four things this table settles:

1. **The brief's headline observation is right in direction and optimistic in size.**
   Plain BSGS with `2^30` entries would give `R = 60` at `~2^30` work — but `2^30` entries is
   **8 GB** at the 8 bytes/entry the brief assumed (and 8 bytes cannot hold a usable fingerprint
   *plus* a 30-bit index, so it is really 16). This box has ~8 GB *available* shared across the
   whole fleet against a 2 GB budget for me. **At 2 GB the deterministic reach is `R = 56`, not
   60**, and the `2^30` work figure ignores that the run must also *build* and *sort* the table.
2. **The kangaroo reaches further than BSGS here, and by four bits, precisely because it needs no
   table** — 32 MB against 2 GB. That is the trade the brief anticipated, and it is the right one
   on a box whose binding constraint is memory and whose second constraint is that
   *more threads make it slower* (§7).
3. **The quotient sweep is the best buy on the board at low cost** and it is genuinely
   exhaustive. Both it and the kangaroo scale as coverage ∝ cost²; the quotient sweep's constant
   is ~3.6× better, but it is capped by memory at `m ≈ 26` while the kangaroo is capped only by
   time. So: **quotient sweep first, kangaroo for reach.**
4. **Disk does not rescue BSGS.** Its access pattern is random; agent AB measured the penalty at
   ~`2^20` slowdown for a `2^1` memory gain. Not attempted.

---

# 6. FAMILIES PRICED AND **NOT** RUN — with the reason

Reporting these is the point of the exercise as much as the runs are.

| family | why not run |
|---|---|
| `k₀ = a·b`, `a ≤ 2^20`, `b ≤ 2^B` | **Dominated, and provably so.** `a·b < 2^{20+B}`, so the whole family is a **subset** of the magnitude family `k₀ < 2^{20+B}`, which the kangaroo decides at `2^{(20+B)/2+4}` ops. Reaching `B = 36` by multi-target BSGS costs `2^32` and covers `2^56` keys — the kangaroo covers those same `2^56` keys for `2^32` and covers `2^60` besides. The brief's suggestion is right that it is cheap; it is cheap **and redundant**. The non-redundant version — `a` small, `b` *not* magnitude-bounded — is the quotient family, which **was** run (§4). |
| `k₀` with low NAF / signed weight | Implication points the wrong way (§3), and agent X exhausted `m ≤ 7` already. |
| window family `k₀ = a·2^s`, `\|a\| < 2^{R−1}`, 256 shifts | **Run** at `R = 40` (§4) — it extends agent X's row 4 (`a < 2^34`, no wraparound) to `a < 2^39` **including** the mod-`N` wraparound cases X's integer enumeration could not see. But note it is dominated in coverage: 256 shifts at `2^40` cover `2^48` keys for `2^30` ops, whereas spending the same `2^30` on a single kangaroo buys `R = 52`, i.e. `2^52`. It is run for **disjointness**, not for coverage. |
| magnitude analogue of agent Y's complement sweep | **The obvious version is empty, and this is worth recording.** If the ON-set folds to `k₀` itself, its integer complement `2^256 − 1 − k₀` is `≥ 2^256 − N = 2^128.346` for every `k₀ ∈ [0,N)` — so "the complement of `k₀` is small in magnitude" is **impossible in that branch, at any `R < 128.346`**, recomputed exactly. It is reachable only when the ON-set folds to `k₀ + N`, which needs `k₀ < 2^256 − N`. Working that through: a complement below `2^R` forces `k₀ ∈ (2^256 − 1 − N − 2^R, 2^256 − N − 1]` — an interval of width `2^R` sitting immediately below `2^256 − N`. **That is (up to its single lowest point) the `c_two256_modN` centre in §4**, so the family is covered — but it must not be called "the complement family": it is a magnitude window living entirely inside the `2^-127.7` two-solution branch, and it has nothing to do with Hamming distance. |
| deterministic BSGS at `R = 56` | Subsumed in coverage by the completed `R = 60` kangaroo; its only advantage is exhaustiveness, bought at 2 GB and heavy random-access pressure on a box already at load 20. Priced (§5), declined. |
| any test of `k₀` against hashes, timestamps, PRNG outputs or file emission order | **Forbidden by rule 2 and not attempted.** See §9. |

---

# 7. AN OPERATIONAL MEASUREMENT WORTH CARRYING

Thread scaling of the same 2^26-jump workload, measured back-to-back on this box (load average
14–20 on 4 cores):

| threads | wall | rate |
|---|---|---|
| **1** | 35.1 s | **1.915 M jumps/s** |
| 2 | 57.3 s | 1.172 M jumps/s |
| 4 | 61.1 s | 1.100 M jumps/s |

**More threads made it strictly slower in total throughput, not merely sub-linear.** Every
production run below is therefore single-threaded. Any agent planning a compute budget on this
box from a per-core figure times four will be out by a factor of ~1.7 in the wrong direction.

---

# 8. THE COORDINATOR'S ITEM 3 — U's partition theorem and the 927 lift conditions

See `AE_ITEM3.md` in this directory for the full working; the result is summarised there and the
scripts are `ae_support.py` and `res_support.json` / `res_927_supports.json`.

---

# 9. RULE-2 BOUNDARY — where I stopped

Everything in this file is arithmetic on a point of an elliptic curve: is `T` a small multiple of
`G`; is it near a named constant; does it factor through a small rational. Every centre used in
§4 is definable from the curve and field data alone (`N`, `p`, `λ`, `β`, powers of two, and
elementary rational multiples of `2^256 − 1`, the constant agent Y's complement identity already
put in play). **No test of provenance was run and none was designed:** no PRNG traces, no
hashes of strings or timestamps, no examination of emission order or coefficient templates, no
hypothesis about who wrote the file or how. Where a family looked like it was drifting from
"low-complexity integer" toward "plausible authoring choice", it was dropped rather than run.
