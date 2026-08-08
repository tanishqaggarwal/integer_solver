# S13 part 3 — THE COMPLETE QUBO SET

Everything below is produced and measured by `s13/qubo_full.py`
(`verify` = correctness, `census` = the block set). Block types are brute-forced
against ground truth, not asserted.

---

## 1. The three primitive blocks (all VERIFIED sound + complete)

Every constraint in the instance is one of three shapes. Each is a real QUBO
(`{(i,j): coeff}`) whose **minimum energy is 0 exactly on the satisfying
assignments** — checked by exhaustive brute force:

```
MUL w=2 q=3 : valid triples   9 / zero-energy   9   sound=Y complete=Y  coupler 7 bits
MUL w=3 q=5 : valid triples  25 / zero-energy  25   sound=Y complete=Y  coupler 7 bits
MUL w=3 q=7 : valid triples  49 / zero-energy  49   sound=Y complete=Y  coupler 7 bits
MUL w=4 q=11: valid triples 121 / zero-energy 121   sound=Y complete=Y  coupler 7 bits
MUL w=4 q=13: valid triples 169 / zero-energy 169   sound=Y complete=Y  coupler 7 bits
LIN w=3 q=5  [1,1]    : sound=Y complete=Y   coupler  9 bits
LIN w=3 q=7  [2,3]    : sound=Y complete=Y   coupler 11 bits
LIN w=4 q=11 [3,5,7]  : sound=Y complete=Y   coupler 13 bits
=> PASS -- every block type is sound and complete
```

**The design decision that makes this work.** The naive encoding
`(x·y − out − k·q)²` is *sound but incomplete*: its column sums go negative while
binary carries cannot, so valid assignments have no zero-energy representation.
I hit this bug and fixed it by writing the identity as `x·y = out + k·q` and
building **two non-negative carry chains** that are then equated bit by bit. Both
sides only add, so every carry is a non-negative binary word. Completeness went
from 2/4 to 5/5, and the max coupler dropped from **5,120 (13 bits) to 80
(7 bits)** — and, importantly, the MUL coupler is **independent of the modulus**.

| primitive | size | coupler |
|---|---|---|
| MUL sub-block (w=16) | **681** binary (excl. shared wires) | 1,280 (**11 bits**) |
| LIN block (w=16) | **362 fixed + 13.8/term** (measured T=2,4,8,16) | ≤ 13 bits |

---

## 2. The complete block set, four tiers

`qubo_full.py census` at `q = 65521`, `w = 16`, target band 1,000–5,000:

| tier | what is unknown | atoms | MULs | total binary | blocks | in band |
|---|---|---|---|---|---|---|
| **[A]** full instance | everything | 42,267 | 253,421 | 189,647,046 | 39,909 | **100%** |
| **[B]** core-anchored | the residual | 5,317 | 5,300 | 5,716,077 | 1,205 | **100%** |
| **[B′]** constant-folded | 328 driver descendants | 433 | 1,010 | 879,087 | 187 | 186/187 |
| **[B″]** branched | `x8731, x9118` only | 229 | 187 | **234,968** /branch | **50** | **50/50** |
| **[C]** p-divisibility core | — | — | — | **79,500** | 32 | ✓ |

Each tier is a *sound* encoding; they differ only in how much verified structure
is assumed fixed. Tier [A] assumes nothing (and is a scale statement, not a
machine target: ~3.8×10⁹ binary over 20 primes). Tier [B″] fixes the 11 already-solved
advice numbers and enumerates the 2 selector bits outside the QUBO (only 4 cases).

**Coupling:** 5,359 shared wires × 16 bits = 85,744 binary of interface against
5,716,077 internal — **1.5%**. Blocks touch each other only through w-bit wire
residues; there is no other interaction.

**Against the stated budget:** the request was blocks of ~1,000–5,000, extendable
to 100,000. Tier [B″] is **50 blocks, max 4,997, all in band**, and the whole
p-divisibility core [C] is **79,500 — a single sub-100k QUBO**.

---

## 3. What the set actually encodes

**Layer A — the circuit, in RNS.** The file is a pure polynomial system over ℤ, so
reduction mod a small prime `q` is a ring homomorphism and every equation has an
exact mod-`q` image. One block per atom; wires are `w`-bit residues shared between
blocks. Disjoint across primes.

**Layer B — the p-divisibility core.** Mod `q`, every p-handle absorbs its check
(`p` is invertible mod `q`), so the RNS layer is *blind* to the real content. That
content is exactly two divisibilities (`core_print.py`):

```
p | x7075·x9118        and        p | x7075·x8731 ,      x7075 = 1 − x2081·x4287
```

emitted as limb/carry chains, verified sound+complete in `qubo_limb.py`.

> **The two layers are complementary and neither is redundant:** RNS handles all
> 39,033 equations cheaply but cannot see divisibility by `p`; the limb chains see
> `p` but only for the core. Together they cover the instance.

---

## 4. Honest limits

- Tier [A] is not a machine target. The practical set is [B″] + [C].
- Every claim here is about **encoding**, not searchability. The blocks are small,
  low-precision and minimally coupled; whether annealing finds their joint ground
  state is untested and not claimed.
- Blocks share wires, so solving them independently is not the same as solving the
  instance — a consensus schedule (block coordinate descent / ADMM over the 1.5%
  interface) is required, and its convergence is likewise untested.
- The deliverable assignment remains **39,026 / 39,033**, unchanged.
