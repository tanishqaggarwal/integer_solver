# S13 — The RNS reformulation: one dense 512-bit QUBO → many disjoint small ones

*Session 13. Responds to the standing question "turn this into disjoint small QUBOs."
Supersedes `REDUCED_PROBLEM.md` §4 where it declares annealing categorically wrong:
§4's **size** and **precision** walls are both consequences of insisting on a single
mod-`p` QUBO, and both dissolve under the reduction below. §4's deeper point — that the
genuine hardness is arithmetic, not combinatorial — survives, and is sharpened here.*

Everything marked **[measured]** was run this session against the file; see
`s13/rns_reduce.py` and the transcript at the bottom.

---

## 0. The one-sentence reformulation

> Because the instance is a **pure straight-line polynomial system over ℤ** (operators are
> only `+ − *` and integer literals — **[measured]** `grep -c '[/%]' EQUATIONS.txt` = 0),
> reduction modulo an integer is a **ring homomorphism**
>
> ```
>        ℤ  ⟶  ℤ/q₁ × ℤ/q₂ × … × ℤ/q_m           (q_j distinct small primes)
> ```
>
> and by CRT an assignment satisfies every equation over ℤ **iff** its residues satisfy
> every equation mod each `q_j`, provided `∏ q_j > 2·max_e |F_e|`. The single dense
> 256-bit-precision optimization problem therefore splits into `m` **independent** problems
> that **share no variables** and whose **coefficients are `< q_j`** — small, and tunable
> independently of `p`.

This is the "another problem that cleanly pulls apart into disjoint small QUBOs" that was
asked for. The disjointness is *exact and structural* (a ring isomorphism onto the product),
not heuristic.

---

## 1. Why the two annealer walls of §4 were artifacts of a single modulus

`REDUCED_PROBLEM.md` §4 rejected annealing on two grounds. Both assumed you must build **one**
QUBO whose objective is `(Σ aᵢbᵢ − t − k·p)²` with 256-bit coefficients.

**§4(a) Size.** "One 256-bit modular multiply ≈ 65,000 binary variables." — True *for a
single mod-`p` multiply*. Under RNS the same multiply is done `m ≈ 20` times, each modulo a
16-bit prime: a 16×16-bit multiply is ≈ 256 partial-product bits, not 65,000. The arithmetic
width per system drops from 256 bits to `⌈log₂ q_j⌉` bits, and the `m` systems are disjoint.
Total logical size is comparable, but it is now a **union of disjoint sparse blocks**, not one
dense blob — exactly the regime hardware wants.

**§4(b) Precision — "the real blocker."** "256-bit coefficients need ≈ 512 bits of coupler
dynamic range; annealers deliver 4–6 effective bits." This is the load-bearing objection, and
it is **the one RNS kills outright.** In the reduced systems every coefficient is a residue
`mod q_j`, so a squared penalty `(… )²` has couplers of magnitude `≤ q_j²`, i.e. width
`≤ 2·log₂ q_j` bits. **[measured]** the failing residual is detected correctly with couplers
as small as:

| modulus `q` | word size | couplers `≤` | catches the 7-defect? |
|---|---|---|---|
| 31  | 5 bits | 10 bits | **7/7**, 0 spurious |
| 61  | 6 bits | 12 bits | **7/7**, 0 spurious |
| 101 | 7 bits | 14 bits | **7/7**, 0 spurious |

Coupler width is now a **free parameter** `2·log₂ q`, traded against the number of systems
`≈ 297/log₂ q`. It is **no longer pinned to 512 bits by `p`.** The precision wall was a
property of the *representation*, not of the *instance*.

**§4(c) "No low-precision residue to peel off."** Correct and unaffected: the 256 message bits
do not control the verification residue (rank 1 of 2). RNS does not claim they do. It attacks a
different axis — the **prime**, not the bit — and that axis is where the separability lives.

---

## 2. Faithfulness, measured on the verified witness

Take the checker-verified `best/new_instance_partial_39026.json` (7 failing over ℤ:
`[12231, 12270, 12350, 14584, 18673, 22044, 29125]`). Reduce every `x_i` mod `q` and
re-evaluate mod `q`. **[measured]** (`s13/rns_reduce.py`):

```
        q  bits   #eqs ≠ 0 (mod q)   caught/7   spurious(false-positive)
        7     3            5            5/7           0
       13     4            6            6/7           0
       31     5            7            7/7           0
       61     6            7            7/7           0
      101     7            7            7/7           0
      251     8            7            7/7           0
     1021    10            7            7/7           0
    65537    17            7            7/7           0
```

Two facts, both important and both honest:

1. **Zero spurious failures at every prime.** A satisfied equation is `0` over ℤ, hence `0`
   mod every `q`. The reduction never *invents* a defect. So mod-`q` satisfaction is a sound
   necessary condition, always.
2. **A single tiny prime is only a probabilistic filter.** `q=7` misses 2 of the 7 (their
   integer residual happens to be `≡0 mod 7`); `q=13` misses 1. This is not a bug — it is the
   CRT redundancy that makes the *ensemble* exact. A prime `q` misses a defect of value `V`
   only if `q | V`; **[measured]** `gcd` of the seven failing values is `1`, so no prime
   misses all seven, and any `q ≳ 31` already catches all seven here. The failing residuals
   are ≈ 2450-bit numbers, so "accidentally divisible by `q`" has probability `≈ 1/q` per
   equation — vanishing across an ensemble.

**Reconstruction schedule.** The genuinely-unknown planted quantities are the thirteen
≈296-bit numbers (the handles are a closed-form `/p` lift, not solved for). To CRT-recover a
296-bit integer you need `∏ q_j > 2^297`, i.e. **[measured]**:

| prime size | # disjoint systems | couplers ≤ | note |
|---|---|---|---|
| 16-bit | **~20** | 32 bits | 3030 primes available — comfortable |
| 12-bit | ~27 | 24 bits | 255 available |
| 8-bit  | ~54 (span 2..256) | 16 bits | only 23 in the top octave; widen the band |
| 5-bit  | ~60 (with prime powers) | 10 bits | annealer-precision-clean |

So: **~20 disjoint systems at 16-bit words with 32-bit couplers, or ~60 disjoint systems at
5-bit words with 10-bit couplers** — pick the point on the curve your hardware's precision
allows. Either way the 512-bit wall is gone.

---

## 3. The second decomposition axis: the circuit DAG (small, sparse, local)

RNS makes the systems disjoint *across primes*. Within one prime the system is still the whole
circuit (≈38,748 residue variables over GF(`q`)) — disjoint but not yet *small*. The circuit
structure supplies the second cut, and it is already mapped in this lab:

- **[prior, s9]** oriented circuit = **31,475 gates** (each an atom `out = f(in₁,in₂,…)` of
  degree ≤ 4, output coefficient ±1) + **7,273 free inputs** + **10,792 checks**.
- Over GF(`q`) each gate is a **tiny constraint on ≤ ~5 residues of `log₂ q` bits each** —
  e.g. `out ≡ in₁·in₂ (mod q)` is a `(out − in₁·in₂ − k·q)²` gadget on `~5·log₂ q + log₂ q`
  binary variables. For `q` at 8 bits that is a **~50-variable QUBO block**.
- The blocks are chained only along **DAG edges** (a residue variable is shared by the O(1)
  gates on its wire). So each per-prime system is a **sparse network of ~50-variable blocks**,
  local coupling, low precision — the native annealer regime — instead of §4's dense
  512-bit-coupler monolith.

Three families of blocks are literally disjoint and literally small, no chaining:

1. **[measured, s9]** the **297 tiny homogeneous components** — each its own small block.
2. **[prior]** the **per-bit load-pin gadgets** — each message bit `b` sits in
   `b·(x_B − HUGE) − s·x_C`; that is a self-contained 1-bit-plus-a-few-residues block.
3. **[prior]** the **residual arithmetic core**, `REDUCED_PROBLEM.md` §3:
   `S = A·u² − w²`, `T = B·u − w·c`, requiring `S ≡ T ≡ 0`. Five unknowns. Over the *small*
   field GF(`q`) this is a **brute-forceable** block (enumerate, or one easy square root —
   Tonelli–Shanks is polynomial for any modulus). The Legendre *existence* obstruction that
   closes the branch mod `p` (`A·c² ≡ B²` needs `A` a QR) is a **single-prime** phenomenon; in
   the RNS ensemble each `q_j` sees its own independent quadratic residue question, and the
   CRT glue — not any one prime — carries the arithmetic content.

---

## 4. What this does and does not buy (honest ledger)

**Buys:**
- **Exact disjointness across ~20–60 primes** — a genuine product decomposition, provable
  (ring isomorphism), not a heuristic partition. **[measured]** validity + faithfulness.
- **Coupler precision decoupled from `p`** — the §4(b) "fatal" wall becomes a tunable. This
  overturns the categorical "annealing is the wrong machine" verdict on precision grounds.
- **Small, sparse, local per-prime blocks** — via the DAG cut — matching hardware topology.
- **A brute-forceable residual core per prime** — the 256-bit square-root miracle needed mod
  `p` is replaced by easy small-field arithmetic plus CRT.

**Does not buy (and §4's surviving truth):**
- RNS does not make the *combined* problem trivial. Each per-prime system must still be
  **solved** (find residues satisfying its checks); the cross-prime agreement that a *single*
  integer solution exists is pushed into the **CRT-lift consistency step**. For a planted,
  essentially-determined instance each prime's solution should be unique and lift cleanly, but
  that is a **[conjecture]**, to be tested (§5), not a theorem.
- It does not manufacture combinatorial control the bits do not have (§4(c) stands).
- Bad primes (those dividing a check's discriminant/resultant, or a defect value) must be
  skipped; the ensemble is over-provisioned to absorb them.

---

## 5. Next experiments (highest EV first)

1. **Per-prime end-to-end solve.** Port the mod-`p` forward evaluator (`s11/gmp1.py`, 0.08 s)
   and the advice-DAG Gauss–Seidel (`s10/advgraph.py`) to a **generic small modulus `q`**.
   For `q = 65521` run: pin the 4 constant numbers, Gauss–Seidel the 9 advice numbers mod `q`,
   forward-evaluate, brute-force the 5-unknown residual core over GF(`q`). Success criterion:
   all 39,033 checks `≡ 0 mod q`.
2. **CRT lift.** Do (1) for ~20 sixteen-bit primes; CRT-combine the thirteen numbers; verify
   the lift with the exact-ℤ `checker.py`. This is the direct route to a **full solve** that
   never touches 256-bit-precision search — the arithmetic is always small-word.
3. **Consistency screen.** If a prime yields multiple residual-core solutions, keep all
   branches and let CRT prune: only the branch agreeing across primes lifts to an integer.
   Measure the branching factor — it bounds the real search cost.
4. **Feed the annealer the honest object.** Emit one per-prime block network (family 1–3 of
   §3) as a QUBO at `q = 251` and hand it to simulated / quantum annealing as a proof of scale:
   ~50-variable blocks, ≤16-bit couplers, sparse — a fair test that §4's precision verdict no
   longer applies.

---

## Artifacts

- `s13/rns_reduce.py` — loads and compiles the system, proves the reduction faithful on the
  verified witness across eight primes, prints the CRT schedule. Runs in ~16 s.
- Validation transcript: satisfied 39,026 equations reduce to `0` mod every prime tested;
  the 7-defect is recovered exactly for every `q ≥ 31`, with **0** spurious failures anywhere.
