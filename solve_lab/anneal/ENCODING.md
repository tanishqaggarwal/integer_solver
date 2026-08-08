# Encoding `EQUATIONS.txt` on a quantum annealer

Everything below was derived and measured in this directory. No step assumes what
the instance "is"; the structure is recovered from the file's own constants and
then verified arithmetically.

---

## 1. What has to be encoded

`EQUATIONS.txt` is 39,031 polynomial equations over `x_0 … x_38747`. All but a
256-bit kernel is solved deterministically by propagation (earlier sessions,
`METHOD_SUMMARY.md`). `reduce.py` and `structure.py` collapse the kernel from the
raw constants, checking each claim:

| step | check | result |
|---|---|---|
| the 256 boolean selectors each gate exactly two ~296-bit constants | `pinrec.json`, dedup by target | 256 pairs |
| a single cubic `y² = x³ + a₂x² + a₄x + a₆ (mod p)`, `p = 2²⁵⁶−2³²−977`, fits **all** 256 pairs | 3-point fit, then verify | 256/256 on the curve |
| the target pair `huge_consts.json` lies on the same cubic | direct | yes |
| depress `x → x − a₂/3` | exact | `y² = x³ + B`, **A = 0** |
| `B/7` is a sixth power mod p | `nthroot_mod` | yes ⇒ F_p-isomorphic to `y² = x³+7` |
| `n·P = O` for all 257 points, `n` = secp256k1 group order | scalar multiply | 257/257, and `n` is prime |
| the 256 points form one doubling chain | `2·P_i == P_{i+1}` | chain length 256, `P_i = 2^i·P_0` |

So, with the selectors written `b_i`:

> **find `b_0 … b_255 ∈ {0,1}` with `Σ b_i·(2^i G) = T` on `E(F_p)`**
> — equivalently `k·G = T` for the 256-bit integer `k = Σ b_i 2^i`.

The group order is **prime**, so there is no Pohlig–Hellman split and no
small-subgroup shortcut: the multi-execution schemes below cannot come from
factoring the group. That was the first thing worth checking, because a smooth
order would have made a minimal annealer encoding trivial.

---

## 2. What "minimal" can mean here

An annealer minimises a quadratic form in binary variables. The 256 answer bits
are an information-theoretic floor. Everything above them is the *verifier*: the
energy function has to be able to tell `k` from every other value, and "is
`k·G = T`?" is not a quadratic function of the bits. So

> minimal encoding = minimal arithmetic circuit for the verifier, compiled to a
> quadratic form with as few ancillas and as little coupler precision as possible.

Six choices do the work, each measured:

**(a) Signed / offset digits — no conditional additions.** Writing
`e_i = 2b_i − 1` (or, windowed, digits `1 … 2^w` with the offset folded into the
target classically) makes every ladder step an *unconditional* addition. No
identity element ever appears, no MUX between "add" and "skip", and at `w = 1`
the addend's x-coordinate is a compile-time constant.

**(b) Comb windows with one-hot digits.** All multiples `2^{wj}·G` are
precomputed classically, so a window's addend is a one-hot *linear* combination
of constants — a table look-up costs no multiplications at all. Window width `w`
trades `⌈256/w⌉` group additions against a `2^w`-entry look-up. Measured optimum
is `w = 8–10` (§4).

**(c) Affine addition with three multiplications.** The circuit inside
`EQUATIONS.txt` uses the fraction-free pair (5 multiplications, plus an
invertibility check to close the `x₂ = x₁` degeneracy — the same "weak division
wire" defect noted in `CIRCUIT_STRUCTURE.md`). Instead:

```
lam·d  = e                 d = x₂−x₁,  e = y₂−y₁
lam·lam = x₃ + x₁ + x₂
lam·(x₁−x₃) = y₃ + y₁
d ≠ 0 (mod p)              two "word ≠ constant" gadgets, ~2·log₂(s) ancillas
```

Three multiplications instead of four (an explicit `d⁻¹` witness) or six. The
non-degeneracy gadget is `hamming(d, c) = 1 + slack`, which is linear — it costs
9 ancillas rather than a whole modular inverse. **Measured saving: 27%**
(`window256.json` vs `window256_neq.json`).

**(d) Never materialise a product.** Each relation is asserted as an exact
integer identity `LHS − RHS − p·q = 0` with an explicit quotient word `q`.
Nothing is ever reduced mod p as a separate step, and no intermediate product
word exists.

**(e) Column balancing, never `(big linear form)²`.** Squaring a form with
`2^{512}` coefficients would need astronomical coupler precision. Instead every
identity is balanced column by column with bounded carries, so all coefficients
stay small.

**(f) Two carry disciplines, a real trade-off.**
`binary` encodes each column's carry as a small binary integer (fewest qubits);
`wallace` compresses columns with 3:2 adders so every penalty is a square of a
form with coefficients in `{±1, ±2}` (lowest dynamic range). At 256 bits:
binary is 2.6× smaller but demands **2²¹–2²³** coupler resolution; wallace holds
at **2⁹**. Neither is inside a real annealer's ~4–5 bits, but wallace is only
about 4 bits over.

Penalties are sums of squares plus AND penalties, so `E ≥ 0` always, and `E = 0`
exactly on solutions. The AND weight is scoped to AND variables' own coupling
load, which alone bought 2^13 → 2^6 at small `s`.

---

## 3. The encoding is faithful — verified, not argued

`demo.py` / `demo_win.py` build the real Hamiltonian for scaled instances and
then **enumerate every candidate scalar**, filling all ancillas by replaying the
construction:

```
[w=2 neq] p=97 m=4 windows=2: 1039 vars, 3725 couplers
    zero-energy k: [9]  true solutions: [9]  -- faithful
```

For every `(w, mode, neq)` combination the zero-energy set is exactly the set of
true discrete logs — no spurious ground state, and the true one always reaches
`E = 0`. The degenerate-division loophole is closed: without the `d ≠ 0` gadget
an assignment with `x₁ = x₂` leaves `x₃, y₃` free and can teleport to the target.

One completeness caveat: the encoding assumes the true `k`'s ladder never hits
`S_i = ±addend` or `O`. At 256 bits that is ~`2^-247`; if it ever mattered, a
different digit offset re-randomises the chain.

---

## 4. Measured cost of the real instance

Marginal cost of one comb window at `s = 256`, built for real (not extrapolated) —
`measure256.py`, `measure_w.py`; a fitted model over `s = 16…64` agrees within 5%.

Run `python3 report.py` for the live table. Summary at the optimum:

| encoding | window `w` | windows | **logical qubits** | couplers | coupler `|J|` range |
|---|---|---|---|---|---|
| `binary` (qubit-minimal) | 9 | 29 | **9.06 × 10⁶** | 3.3 × 10⁸ | 2²¹ |
| `wallace` (precision-minimal) | 8 | 32 | **2.66 × 10⁷** | 1.2 × 10⁸ | 2⁹ |

(`w` swept 1…11 for `binary`, 1…10 for `wallace`. Both curves are flat near the
optimum — binary sits at 9.1–10.4 × 10⁶ for `w = 7…11`, wallace at 2.7–3.1 × 10⁷ for
`w = 7…9` — and turn up as the `2^w`-entry look-up starts to dominate the four
multiplications it saves.)

Smallest indivisible piece: **one 256×256 modular multiplication ≈ 8 × 10⁴
(binary) to 2.2 × 10⁵ (wallace) qubits.** Nothing smaller is a self-contained
subproblem.

Against hardware (D-Wave Advantage 5,760 qubits / Advantage2 4,400, degree 15–20):

* **1,600× – 6,000× short on qubits**, before minor-embedding overhead. The
  measured average degree is ~9, which is under the hardware degree, so chains
  would be short — embedding is not the bottleneck, raw count is.
* **4 to 16 bits short on coupler precision.**
* Largest instance of this exact shape that fits today: a **~13-bit curve**
  (`extrapolate.py`). The instance is 256-bit.

---

## 5. Multiple executions: what each scheme actually buys

| scheme | per-run qubits | runs | sound? | verdict |
|---|---|---|---|---|
| **D1 interval split** — enumerate the top `256−μ` bits classically (a classical EC subtraction per run), anneal the low `μ` | `⌈μ/w⌉ ×` window | `2^{256−μ}` | yes, exact | the only exact decomposition. To beat Pollard rho the anneal must carry `μ ≥ 128`, i.e. **≥ 16 windows ≈ 4.7 × 10⁶ qubits**; below that the classical outer loop *is* the whole search |
| **D2 ladder split / meet-in-the-middle** — anneal the two halves against a shared 512-bit accumulator boundary | half | 2 | only if the boundary is enumerated | the boundary collision is exactly BSGS: `2^128` time *and* memory |
| **D3 hybrid LNS / qbsolv** — clamp all but a sub-block, anneal, iterate | any size you like | many | no | measured to fail (§6) |
| **D4 modulus relaxation** — assert the identities mod a 32-bit prime instead of over ℤ, intersect across runs | ~`3.4 × 32²` per multiply | many | no | shrinks a modmul ~64× but admits ~`2^{256}/m²` spurious solutions; the runs cannot be intersected classically |

**The floor.** Multiple executions cannot push the per-execution footprint below
one modular multiplication (~10⁵ qubits) without giving up soundness, and cannot
push the *useful* footprint below ~16 windows (~5 × 10⁶ qubits) without the
classical outer loop doing all the work. That is the honest answer to "minimally,
with multiple internal executions".

---

## 6. Does the landscape anneal at all?

Qubit count is not the only obstacle. `unit_probe.py` takes the **atom** of the
encoding — a single modular multiplication — and asks simulated annealing to
settle it (8,000 sweeps × 8 restarts):

```
  s    p    vars  coupl |  a,b clamped (ancillas only) |  free search
  6   41     153    543 |  3/8 reach E=0               |  0/8, best E=1
  8  163     255    970 |  2/8 reach E=0               |  0/8, best E=2
 10  641     359   1373 |  0/8, best E=3               |  0/8, best E=5
 12 2539     627   2571 |  0/8, best E=9               |  0/8, best E=10
```

`sa_probe.py` does the same to a whole 4-bit ladder (2,726 variables) with the
**answer digits clamped to the known solution** — the easiest possible task on
this Hamiltonian, a unique forced completion with no search over `k` at all:

```
  2,000 sweeps: 0/6 reach E=0 (best E=74)
 20,000 sweeps: 0/6 reach E=0 (best E=50)
100,000 sweeps: 0/6 reach E=0 (best E=40)
```

A 50× budget increase moves the best energy 74 → 40. The annealer cannot even
*fill in ancillas whose values are uniquely determined*, let alone search. This
is the expected shape for a group-hard problem: flipping any subset of answer
bits sends the accumulator to a pseudorandom point, so the energy carries no
gradient toward `k` — a golf course, not a funnel. It also explains why D3
(hybrid decomposition) cannot work: block-coordinate descent needs a landscape
where local improvement means global progress, and there is none here.

---

## 7. Bottom line

The minimal faithful annealer encoding of `EQUATIONS.txt` is
**≈ 9.1 × 10⁶ logical qubits** (qubit-minimal, 2²¹ coupler precision) or
**≈ 2.7 × 10⁷** (precision-minimal, 2⁹). Both are verified exact: ground states
are precisely the solutions of the original 39,031 equations.

That is ~10³–10⁴× more qubits than exists, 4–16 bits more coupler precision than
exists, and — the part extra hardware would not fix — a landscape on which
annealing measurably does no better than guessing. Multiple executions move the
per-run footprint but not the product: the exact decomposition (D1) is a
classical exhaustive search wrapped around the QPU, and the schemes that shrink
the QPU's share below one modular multiplication all give up exactness.

---

## 8. Files

| file | role |
|---|---|
| `reduce.py` | 39,031 equations → curve + 256 gated points + target (`core.json`) |
| `structure.py` | verifies the curve, the prime order, the doubling chain; states the core problem |
| `qubo.py` | QUBO compiler: AND cache, column balancing, both carry disciplines, witness replay |
| `ladder.py` | the encoding: signed-digit ladder (`build`) and windowed comb (`build_win`) |
| `demo.py`, `demo_win.py` | exhaustive faithfulness proof on scaled instances |
| `resources.py`, `measure256.py`, `measure_w.py` | marginal window cost, measured at `s = 16 … 256` |
| `report.py`, `extrapolate.py` | resource tables, optimal `w`, hardware comparison |
| `sa.py`, `sa_probe.py`, `unit_probe.py` | the annealability experiments |
| `core.json`, `window256*.json`, `unit_probe.json` | raw data behind every number above |

```bash
python3 reduce.py && python3 structure.py     # the reduction, re-derived
python3 demo_win.py                           # faithfulness, exhaustive
python3 report.py                             # the resource table
python3 unit_probe.py                         # can it anneal?
```
