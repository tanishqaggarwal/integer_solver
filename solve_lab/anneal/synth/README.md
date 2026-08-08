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

## The two open levers (agents in `synth/mincost/`, `synth/solver/`)

- **mincost**: compose the measured wins (Toom-3 modmul, pseudo-Mersenne
  reduction, MUX/AND-tree one-hot, signed digits, Montgomery x-only ladder) to
  shrink one window ~3×, which raises every μ above.
- **solver**: whether a real annealer *converges* on this landscape at all. Prior
  finding: plain SA cannot pass ~8-bit modular multiplies. If the true cap on
  bits-per-run is solver reach rather than qubit count, that — not the table
  above — is the binding constraint. This track charts it honestly.

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
