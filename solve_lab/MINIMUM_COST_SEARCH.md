# Minimum-cost search to solve the instance

Coordinator synthesis over the whole campaign. Every input is a measured result recorded in
`FLEET.md` or an agent's RESUME file; the arithmetic below was computed, not estimated.

## 1. What the problem reduces to

Established by **four independent parses sharing no code** (agents F, K, L, Q, corroborated by P,
R, T, U, W):

> The instance is satisfiable **iff** some subset `S ⊆ {0..255}` of leaf selectors satisfies
> `k·G = T`, where `k = Σ_{i∈S} 2^i`.

Supporting facts, each verified at instance level:

- **The lift is not an obstacle.** Agent T closed the 927 integer conditions at
  `|S| = 2, 3, 5, 6, 7, 8, 17`, each dumped and checked. So **finding `S` yields 39,033/39,033**;
  the remaining work after `k` is minutes.
- **The classification closes at exactly two families** (agent W): the Jacobian
  `d(N1,N2)/d(i5,i6)` has `det = A³`, so `A ≡ 0` forces `B ≡ 0` (degeneracy, output free) and
  `A ≢ 0` determines the output uniquely (chord). **No third case; the mixed case is impossible.**
  *Scope: this is at atom level. Equation-level cancellation is a larger set the theorem does not
  cover.*
- **The degeneracy route is closed** (agent U, exhaustive): all 510 proper slot supports have masked
  value strictly below `N`, largest ≈ 0.7987·N, so no two gadget inputs can coincide by
  configuration.
- **The solution is essentially unique.** `2²⁵⁶ − N ≈ 2¹²⁸`, so there are one or two valid `S`, and
  two only with probability ≈ 2⁻¹²⁸. No multi-target advantage exists.

## 2. The curve is secp256k1 — verified bit-for-bit

    p == 2^256 - 2^32 - 977                      : True
    N == 0xFFFF...FFFEBAAEDCE6AF48A03BBFD25E8CD0364141 : True

This is why every structural shortcut checked in this campaign failed, and it is not a coincidence
of this instance: prime order, **not** anomalous (`N ≠ p`), embedding degree past 24 (verified to
k ≤ 24), and the only endomorphism is the standard order-3 `λ`. Fifteen years of public
cryptanalysis have produced no better-than-generic attack on it.

## 3. The minimum-cost known attack

**Parallel Pollard rho with distinguished points, negation map, and the GLV endomorphism.**

| variant | cost |
|---|---|
| plain rho, `√(πN/4)` | 2^127.8 |
| + negation map (`√2`) | 2^127.3 |
| + GLV endomorphism (`√3`) | 2^127.0 |
| **+ both (`√6`)** | **2^126.5 group operations** |

Engineering specification, so this is a real answer rather than a shrug:

- **Iteration**: r-adding walk, r ≈ 20 partitions, walking on equivalence classes under the
  order-6 automorphism group `{±1, ±λ, ±λ²}`.
- **Fruitless cycles**: the negation map's known pathology. Detect by cycle-length check and escape
  by a deterministic perturbation; without this the `√2` is lost.
- **Distinguished points**: `x`-coordinate with `d` leading zero bits. Storage ≈ `2^126.5 / 2^d`
  points; per-worker memory `O(1)`.
- **Parallelism**: linear speedup in workers (van Oorschot–Wiener), no communication except DP
  reports.
- **Termination**: a DP collision yields `k` directly; then run agent T's `t_close2wj.py` lift and
  verify with `checker.py`.

### What it costs in the real world

| machine | wall time for 2^126.5 |
|---|---|
| one CPU core (10⁶ ops/s) | 3.9×10²⁴ years |
| one GPU (10⁹ ops/s) | 3.9×10²¹ years |
| one million GPUs (10¹⁵ ops/s) | 3.9×10¹⁵ years |
| **Bitcoin-network scale (6×10¹⁸ ops/s)** | **6.5×10¹¹ years ≈ 47× the age of the universe** |

**This is not a resource question.** No budget reaches it.

## 4. The only thing that beats rho: bounded-complexity `k`

If `k` is drawn from a small structured class, meet-in-the-middle over that class beats `2^126.5`.

**Unsigned Hamming weight** — `k` has `w` one-bits. MITM cost ≈ `C(256, w/2)`:

| bound | cost | vs rho |
|---|---|---|
| w ≤ 10 | 2^33.0 | beats |
| w ≤ 20 | 2^58.0 | beats |
| w ≤ 30 | 2^79.1 | beats |
| w ≤ 40 | 2^97.8 | beats |
| w ≤ 50 | 2^114.6 | beats |
| **w ≤ 56** | **2^123.8** | **crossover** |
| w ≤ 60 | 2^129.7 | loses |

**Signed-digit weight** — `k = Σ ε_j·2^{e_j}`, `ε ∈ {±1}`, `m` terms. **Strictly more general**: it
contains every low-Hamming-weight `k`, *and* every low-**run-length** `k` (a run of ones is
`2^{a+1} − 2^b`, two signed terms), *and* anything from a short addition-subtraction chain. Cost
≈ `C(256, m/2)·2^{m/2}`:

| bound | cost |
|---|---|
| m ≤ 8 | 2^31.4 |
| m ≤ 10 | 2^38.0 |
| **m ≤ 12** | **2^44.4** |
| m ≤ 14 | 2^50.6 |
| m ≤ 16 | 2^56.5 |

**Signed-digit `m ≤ 10` costs about what unsigned `w ≤ 12` costs and covers a strictly larger,
structurally different set.** Nobody in this campaign has searched it. Agent X has been tasked with
it after its plain-weight sweep.

### Budget → reachable bound

| budget | unsigned weight | signed-digit weight |
|---|---|---|
| 2^47 (this box, ~8 h) | w ≤ 14 | m ≤ 12 |
| 2^60 | w ≤ 20 | m ≤ 16 |
| 2^70 | w ≤ 24 | m ≤ 20 |
| 2^80 | w ≤ 30 | m ≤ 24 |
| 2^90 | w ≤ 34 | m ≤ 28 |
| 2^126.5 | — | **full solve by rho** |

## 5. Already excluded (all mod p; see the standing note on Q's withdrawal and restoration)

`w ≤ 6` exhausted; `w ≤ 7` partially; `k` not in the bottom or top `2^44`; ON-bits not confined to
any 34-bit window; `k` not a small multiple of a ladder point (`m ≤ 10⁷`); not in the endomorphism
orbit; no `k = a + bλ` with `|a|,|b| < 2^21`.

## 6. The bottom line

**If the instance was generated by choosing a random 256-bit `k` and publishing `T = k·G`, it is not
solvable by any known method, and no search this campaign or any successor can run will change
that.** The reduction is complete, the classification is closed, and the residual is a generic
ECDLP on the most-studied curve in existence.

**The minimum-cost search that can actually be run** is therefore bounded-complexity MITM —
signed-digit weight first, since it dominates unsigned at equal cost. It wins **only** if `k` lies
in a low-complexity class, which requires a prior the mathematics cannot supply.

**The one place such a prior could come from is how the instance was constructed** — the
leaf-to-exponent assignment, the ordering of the 512 pin constants, anything reflecting a generator
rather than the mathematics. **That direction was closed by user instruction at the start of the
campaign and the fleet has respected it throughout.** It is the only remaining lever, and it is the
user's to open.
