# The obstruction, reduced to a single conserved functional

This session drove the `EQUATIONS.txt` wall down to its irreducible form. Best verified
remains **39,022 / 39,033** (`best_agentA_39022.json`). The residual is now completely
characterized as a **1-dimensional conserved obstruction**.

## The reduction chain (all verified this session)

1. **The 11 fails are a 2-bit MUX + p-granular slacks.** The only nonzero atoms at agentA are
   - `20862 = x_7068 − x_2099 − 7376877·x_642` (with `x_642 = p·x_17325`)
   - `20864 = x_4432 − x_19964 − x_28730` (with `x_28730 = p·x_9413`)
   so the real conditions are `x_7068 ≡ x_2099 (mod p)` and `x_4432 ≡ x_19964 (mod p)`.
   `x_2099`, `x_19964` are a MUX in bits `a=x_2081`, `b=x_4287`:
   `x_2099 = b(1−a)x_31861 + a(1−b)x_6418 + ab·x_9118`, likewise `x_19964`.
   All four channels were scanned (fails after closing G1/G2): **(1,0)=16**, (0,1)=41,
   (1,1)=46, (0,0)=118. Channel (1,0) (agentA) is optimal.

2. **Closing G1/G2 is free but relocates the residual to 16 "leaf-ripple" equations:**
   `[697,1985,5225,10815,16048,17784,17801,22402,23667,24721,27124,28737,29638,29959,35935,37431]`.
   Setting the free leaves `x_7068=x_2099`, `x_4432=x_19964` (slacks 0) fixes the 11 and
   breaks exactly these 16 (the other half of the verifier).

3. **The core is disconnected from the verifier.** Perturbing every core knob
   (`x_14853, x_12186, x_16742, x_22162, x_30213`) moves **0/16** leaf-ripple residuals while
   breaking S, T and all three loads. So agentA's *degenerate* core (`x_29322=x_3558=0`) vs a
   non-degenerate core is **irrelevant** to the verifier — refuted as an escape.

4. **The 16 leaf-ripple equations depend on only 16 free inputs:**
   `[96, 332, 2081*, 2964, 4432, 7068, 8557, 9280, 10528, 13195*, 18027, 24548, 27711, 27863, 29162, 31687]`
   (`*` = selector bit). Their exact GF(p) Jacobian (w.r.t. **all** 8583 free inputs) has
   **rank 12**; `rank[J|R] − rank[J] = 1`.

5. **The obstruction is a single conserved functional.** There is exactly one linear
   combination `c` of the 16 residuals with `c·(∂R/∂u) = 0` for **every** free input `u`, yet
   `c·R ≠ 0`. It is invariant under all continuous moves. The only two discrete knobs in its
   support are `x_2081` and `x_13195`; all four of their on/off combinations leave the count
   `≥ 16` (flipping `x_13195` breaks 40 more on raw agentA). So the conserved value cannot be
   moved by any continuous step or any local bit flip.

## What this means

Within agentA's forward-construction branch the residual is **immovable** — a topological
invariant fixed by the setter's message constants. Crossing it requires a *global*
reconfiguration of the 256-bit message/bit-pattern (the setter's witness), not any local
move. No polynomial weakness in that global step has been found: the pin CONSTs are
cryptographically random (512 distinct residues, gcd 1), each message word appears in exactly
one pin (no duplication shortcut), and the auxiliary modulus 6672769 is prime and coprime to
p (no CRT shortcut).

The remaining live attack surface is the global message-pattern search — pursued by the
cube/higher-order-differential, non-degenerate-core, and algebraic-elimination agents.

## Tooling added this session
`local_tangent.py` (69-col tangent consistency), `closure_test.py` (coupling does not
localize: 69→839→5851), `expanded_tangent.py` (inconsistent at every finite truncation),
`conserved.py`/`conserved2.py` (the 1-dim conserved obstruction extraction).
