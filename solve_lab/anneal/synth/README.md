# Synthetic ECDLP-on-an-annealer: planted keys, end-to-end, no live target

Everything here uses **synthetic curves with a key generated and known locally**
(`synth/gen.py`). Every instance is solved and verified against its planted key.
No live public key is involved.

## The objective and the one identity that governs it

Minimize the **number of annealer runs** to recover a `b`-bit key. The only sound
multi-run scheme (proved in `multirun/`) is the interval split: the classical
outer loop fixes the top `b − μ` bits, subtracts them from the target, and hands
the annealer a QUBO for the low `μ` bits. Hence

> **runs = 2^(b − μ),  μ = scalar bits resolved per run.**

`μ` is capped by (qubits per QUBO) ≤ (machine budget) and by coupler precision.
So **reducing runs ≡ raising μ ≡ shrinking one comb window** — the run count and
the encoding size are the same lever. Every qubit saved per window halves the
runs per additional bit.

## What is demonstrated (measured, this directory)

**1. End-to-end recovery works.** `synth/demo.py` plants a key, builds the comb
QUBO, solves each sub-instance, recovers the scalar, and checks `k·G = T`:

```
bits  mu/run  runs-to-hit  worst 2^(b-mu)  recovered=planted
  16    16         1           2^0              YES
  24    16         1           2^8              YES
  32    16         1           2^16             YES
  36    16         1           2^20             YES
```

**2. The scheme is sound.** Scanning *all* `2^(b−μ)` prefixes, exactly one yields
the planted key (the "two hits" at some sizes are `k` and `k+n`, identical mod n):

```
bits=16 mu=8 : 256 prefixes, 1 success, recovered = planted
bits=20 mu=10: 1024 prefixes, 1 success, recovered = planted
```

**3. The runs table, on measured encoding costs** (`synth/runs_table.py`,
`synth/window_grid.json`): each comb window does **full s-bit field arithmetic**
(s = key bit-size) no matter how few scalar bits it resolves, so one window costs
5.2k qubits at 32-bit, 19k at 64-bit, 275k at 256-bit. Consequently:

| key size | Zephyr 4.4k / Pegasus 5.8k (5-bit J) | Toshiba SB 1e5 (full, fp) | idealized 1e6 |
|---|---|---|---|
| 32-bit | no fit → 2³² | **1 run** | 1 run |
| 64-bit | no fit | 2⁴⁰ | 1 run |
| 128-bit | no fit | 2¹²⁰ | 2³² |
| 256-bit | no fit | no fit | 2²²⁹ |

**Two walls, in order.** On D-Wave the *coupler-precision* wall (≈5 usable bits vs
the encoding's 2¹⁵–2¹⁹) bites before the qubit wall — no window fits at any size
here. Full-precision, fully-connected Ising machines (Fujitsu DA, Toshiba SB)
clear that wall, and then the *qubit* wall sets μ. Neither machine class brings
256-bit into range: even a hypothetical 10⁶-qubit machine needs 2²²⁹ runs.

## The verdict: two walls, and a run count that will not drop

`synth/solver/` settled the decisive question. The encoded landscape has **no
gradient**: fixing all ancillas, a candidate one digit from the key has the same
energy as one that differs in every digit (correlation 0.06). Parallel tempering
and simulated bifurcation — the barrier-crossing solvers — do no better than SA,
and nothing reaches the ground state past ~6-bit modular multiplies. The solution
is a needle in a flat haystack.

That changes the run-count arithmetic. The interval split gives **outer runs =
2^(b-mu)**, and shrinking the encoding raises mu, cutting outer runs — real, and
the encoding is now ~3x smaller (squeeze/) at the arithmetic floor. But with no
gradient a single anneal finds a sub-instance's needle only by chance, so
**anneals per sub-instance ~ 2^mu**, and

> total anneals = 2^(b-mu) . 2^mu = 2^b,  invariant in mu.

Shrinking the window trades outer runs for inner difficulty at a fixed product.
The annealer gives **no speedup over classical brute force**. So there are two
independent walls, either one fatal:

1. **Size** (`squeeze/FINDINGS.md`): one 256-bit modmul needs >=4x a 4,400-qubit
   machine in partial-product ancillas alone; the full ladder is 1.16e7 physical.
2. **Landscape** (`synth/solver/FINDINGS.md`): even where a sub-instance fits, it
   is a gradient-free needle search costing ~2^mu anneals.

The run count cannot be pushed below the classical search bound. That is the
honest, measured answer.

## The encoding lever (still worth its floor: `synth/mincost/`)

Completing the composition still gives the tightest per-window number (it lowers
OUTER runs and decides what *fits*), even though it cannot beat the landscape:

- **mincost**: compose the measured wins (Toom-3 modmul, pseudo-Mersenne
  reduction, MUX/AND-tree one-hot, signed digits, Montgomery x-only ladder) to
  shrink one window ~3×, which raises every μ above.
- **solver** (done): charted above — no gradient, no speedup. See
  `synth/solver/FINDINGS.md`.

## Files

`gen.py` (planted-key generator, prime order via BSGS), `build_lib.py` (cached
instance library), `solve.py` (interval-split recovery), `demo.py` (the
end-to-end run), `runs_table.py` + `window_grid.json` (measured runs scaling),
`measure_grid.py` (per-window cost vs field size).

```bash
python3 -m synth.build_lib 8 12 16 20 24 28 32   # cache instances
python3 -m synth.demo                             # recover planted keys, count runs
python3 synth/runs_table.py                       # the measured runs-vs-machine table
```
