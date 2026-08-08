# S13 part 2 — Core reduction + minimally-coupled QUBO blocks

*Goal: reduce the instance to a small core, then decompose that core into
minimally-coupled QUBOs of ~1,000–5,000 binary variables each, total ≤ 100,000.*
All numbers **measured** this session (`s13/core_reduce.py`, `core_cone.py`,
`core_print.py`, `core_extend.py`, `qubo_limb.py`).

---

## 1. The core, decoded exactly

Reduction chain, each step measured:

| step | size |
|---|---|
| raw file | 39,033 equations / 38,748 variables |
| atoms | 42,267 |
| circuit orientation (s9) | 31,475 **gates** (evaluated, not unknowns) + 7,273 free inputs + 10,792 checks |
| free inputs carrying information | **13 advice numbers** (~296-bit) + 2 set message bits |
| **genuine unknowns** | **3,330 binary** (13 × 256 bits + 2) |
| residual at the 39,026 witness | **2 nonzero checks** |
| **ancestor cone of those 2** | **17 variables, 13 atoms, 6 free inputs, 2 wide multiplies** |

The cone printed symbolically (`core_print.py`), with `x26064 = p = 2²⁵⁶−2³²−977`
propagated along the identity wire (`x35638 = x18822 = x22665 = x28961 = p`), and
`x1329`, `x10903` free handles:

```
a35759 = 5113045·x7075·x9118 − x29854 ,   x29854 = x1329·p
a35760 = −x10903·x28961 + x31864     ,   x31864 = −x7075·x8731
x7075  = 1 − x2081·x4287                  (x2081, x4287 boolean)
```

so the **entire residual of the instance is two divisibilities**:

> ```
> p | x7075 · x9118        and        p | x7075 · x8731
> ```

with `x9118`, `x8731` **free inputs**. Two doors: drive the selector `x7075` to 0
(set `x2081 = x4287 = 1`), or make both free inputs multiples of `p`.

**Why it is nonetheless hard — measured collateral** (`core_extend.py`): the four
drivers `{x2081, x4287, x8731, x9118}` have 328 descendant variables touching 433
atoms; **109 checks are at risk, 28 self-absorb via a private p-handle, and 81 are
hard constraints.** The cone of those is 6,237 variables / 53 wide multiplies. So
the core is small; its *neighbourhood* is what previous sessions kept paying for.

---

## 2. The decomposition: limb/carry chains, not one dense block

The naive encoding of `a·X − b − k·p = 0` is the single penalty
`(a·X − b − k·p)²` — one dense block with **512 bits of coupler dynamic range**,
the wall `REDUCED_PROBLEM.md` §4(b) called fatal.

Instead evaluate the same integer identity **column by column in radix 2^L**. At
column `j`:

```
s_j = (a·X)_j − b_j − (k·p)_j + c_{j−1}
require  s_j mod 2^L == 0        (result limb vanishes)
c_j = s_j >> L                   (carry to the next column)
```

Every coefficient inside a column is `< 2^L`, so **couplers are ~2L bits**, and
column `j` touches column `j+1` **only through the carry word `c_j`**. The result
is a linear chain

```
[col 0] --c₀-- [col 1] --c₁-- … --c₁₄-- [col 15]
```

whose **separator between any prefix and suffix is a single carry word** — minimal
coupling in the strict graph sense; treewidth is the carry width, not the operand
width.

**Correctness is verified, not asserted** (`qubo_limb.py verify`): the checker
walks the actual carry chain (not `abs(a·X−b−k·m)`), and over six random moduli it
confirms **sound** (every zero-energy state satisfies the congruence) and
**complete** (every solution is reached) — 6/6 PASS.

---

## 3. Measured block sizes — this is the answer to the request

`qubo_limb.py size`, 256-bit operands:

**Unknown × unknown multiply** (the general collateral check):

| limb width | blocks | block size min/mean/max | **total binary** | carry coupling | in 1k–5k |
|---|---|---|---|---|---|
| **L = 16** | **16** | 322 / 2,484 / 4,646 | **39,750** | **18–21 bits** | **13/16** |
| L = 32 | 8 | 1,154 / 4,963 / 8,773 | 39,709 | 34–36 bits | 4/8 |

**Known-coefficient congruence** (`X ≡ 0 mod p` — the core condition itself, which
is *linear* in the bits of `X`, so no multiplier array at all):

| limb width | blocks | block size | total binary | carry coupling |
|---|---|---|---|---|
| L = 16 | 16 | 50 | **800** | 18 bits |

**Budget check against the 100,000 cap:**

- one linear 256-bit congruence = **800 binary** → 125 of them fit;
- one full 256×256 multiply = **39,750 binary** → **two fit** (79,500 < 100,000).

> **L = 16 is the operating point**: 16 blocks per multiply, mean **2,484**
> variables per block, 13 of 16 inside the requested 1,000–5,000 band, coupled by
> **18–21-bit carry words**. Under the 100,000 cap this buys two full 256-bit
> modular multiplies, or one multiply plus ~75 linear congruences.

---

## 4. How the pieces compose

- The **two core conditions** are linear congruences → **800 binary each**, trivial.
- Each **collateral check** that genuinely multiplies two unknowns costs one
  39,750-binary chain of 16 blocks.
- Blocks within a chain couple only by carries (18–21 bits). Chains couple to each
  other only through the **shared driver residues** (`x9118`, `x8731`, and the two
  selector bits) — a handful of 256-bit words, i.e. a narrow interface, so a
  block-coordinate-descent / ADMM schedule over the chain graph is well posed.
- Contrast the whole-instance METIS partition (`decompose_metis.py`, retained for
  reference): 460 blocks, mean 2,998 binary, **25% interface fraction** and a hub
  (`x24453`) touching 286 blocks. The **core-first** route is far better coupled
  because the core is 17 variables, not 38,748.

---

## 5. Honest scope

- The 100,000-variable budget covers the **core conditions plus about two wide
  collateral multiplies** — not all 81 hard checks at once (53 wide multiplies
  ≈ 14M binary if encoded simultaneously). The intended use is
  **core + selected collateral per round**, iterating, not one monolithic QUBO.
- Nothing here claims the instance is solved, or that annealing will find the
  ground state. What is established: the **encoding obstacles** of
  `REDUCED_PROBLEM.md` §4 — dense 512-bit couplers and 10⁵–10⁶-variable monoliths —
  are **removed** by the core reduction plus limb/carry decomposition. Whether the
  resulting landscape is searchable is a separate, open question.
- The two doors (`x7075 → 0` vs `p | x9118, x8731`) were both previously explored
  and both relocate cost; the contribution here is that they are now **exact, tiny,
  and encodable**, with the collateral explicitly priced at 81 hard checks.

## Artifacts

| file | role |
|---|---|
| `s13/core_reduce.py` | reduction chain R0→R5; genuine unknowns = 3,330 binary |
| `s13/core_cone.py` | residual cone: 17 vars / 13 atoms / 2 wide multiplies |
| `s13/core_print.py` | symbolic dump of the core — the two divisibilities |
| `s13/core_extend.py` | collateral: 109 at-risk checks, 81 hard, 53 wide multiplies |
| `s13/qubo_limb.py` | limb/carry chain emitter + **verified** sound/complete; sizing |
| `s13/decompose_metis.py` | whole-instance METIS reference partition (460 blocks, 25% cut) |
