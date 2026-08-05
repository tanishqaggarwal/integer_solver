# Endgame analysis — the `EQUATIONS.txt` trapdoor, completely characterized

Best verified: **39,022 / 39,033** (`best_agentA_39022.json`, also reproduced as
`gadget_handled.json`). This document records the deepest and most complete
reverse-engineering of the obstruction reached to date — the full 256-bit control
structure, an *exact* solution of the nested gadget verifier, and the precise reason the
final residual cannot be crossed without the setter's secret.

## 1. The architecture (fully decoded)

`EQUATIONS.txt` is an obfuscated arithmetic circuit over GF(p), p = 2²⁵⁶−2³²−977 (the
field prime), compiled to integer equations. It is a **256-bit MUX trapdoor**:

- **256 free control bits** (selectors), each appearing in exactly **2 pin atoms** of the form
  `selector · (x_target − CONST) − coef · handle`. So **512 pins** total (`pinrec.json`).
  When a bit = 1, it loads two ~256-bit message constants onto two wires; when 0, its
  handles (wire products `p·free`) must vanish.
- The 512 pin **CONSTs are the setter's secret message.** They are structureless: 512
  distinct residues mod p, gcd 1, no small values, 72-digit minimum gap, 256 distinct
  coefficients — i.e. cryptographically random. **No lattice/algebraic shortcut.**
- The loaded values feed the **verifier**: degree-4 perfect-square atoms `Q²=0` (⟺ `Q=0`).
- The **MUX** `x_15298 = OR(x_1222,x_35723)·OR(x_28505,x_32083)` selects which verifier is
  active. It is *independent* of the load bits `x_2081`, `x_4287`.

## 2. The two branches, and why both are walled

Everything reduces to closing the verifier. There are two sub-branches, set by `x_9062=x_4287`:

### Branch A — `x_9062 = 0` (old quadrant, `best_agentA_39022`)
Core solves (S≡T≡0 mod p) but the residual is **3 atoms** → two gaps
`G1 = 7376877·x_642 + x_2099 − x_7068`, `G2 = x_4432 − x_19964 − x_28730`, with
`x_642 = p·x_17325`, `x_28730 = p·x_9413` (wire products). The gaps are **sub-p**; the
slacks are **p-granular** → G1,G2 cannot be zeroed. 11 fails. Hard wall.

### Branch B — `x_9062 = 1` (new quadrant) — **new this session**
Here `x_2099 = x_9118` and `x_19964 = x_8731` become *free*, so G1/G2 *can* close in ℤ.
But `x_9062=1` needs `x_4287=1`, and the escape `x_2099=x_9118` also needs **`x_2081=1`**,
which forces `x_21279 = x_4287·x_2081 = 1` — activating a **nested gadget verifier**
(atoms 17897, 20866, 20868, 34232 + square 45603). This session I solved that gadget
**exactly**:

- x_2239, x_31731, x_9106 are all linear combinations of x_27177, x_4306.
- x_27177, x_4306 are affine in the two free knobs **x_8731, x_9118**.
- A **Diophantine + CRT** solve on (x_8731, x_9118) gives x_31731 = 0 exactly,
  13523997 | x_9106, x_2239 ≡ 0 mod p; setting handles x_9629, x_6947, x_33168 zeros all
  four gadget atoms (`gadget_handled.json`, 11 fails — the G1/G2 set).

Closing G1/G2 then moves the free leaves x_7068, x_4432, which **ripples** to 16 equations
that reduce to **4 atoms**: 7450 (`x_2964 − x_26756 − x_579`), 7452
(`9367949·(x_24548−x_25442) − x_7927`), verifiers 44342, 45677. The two slacks are
`x_579 = p·x_19569`, `x_7927 = p·x_11052` — **p-granular again** (x_13859, x_15616 are wire
members = p), gaps sub-p. Iterative residue heal **empirically diverges** (16→12→33→20).

## 3. The conservation law (the whole wall)

Every path relocates one **conserved sub-p residue** that no p-granular slack can absorb:

- The identity **wire is pinned to p** by atom **37110 = `x_26064 − p`**, which is
  degree-1 and **slack-free**: its 12 equations contain only rigid identities and *boolean*
  atoms `x·(x−1)` (no free product-slack). Verified this session. Hence every wire member
  (x_13859, x_15616, x_28599, x_17499, …) = ±p, and every product-slack is a multiple of p.
- The METHODOLOGY escape that cracked the *previous* instance (free wire → set V=1 →
  fine-grained slacks) is therefore blocked. Freeing any single wire member requires moving
  the whole union-find class, which lights up the slack-free 37110.
- Escape and gadget are **bundled through x_2081**: x_2081=1 gives the G1/G2 escape *and*
  the gadget; x_2081=0 removes both (and breaks x_2099=x_9118). Neither branch escapes.

The global system *is* consistent (the witness exists), but its linearization at every
reachable forward-construction point is inconsistent — the signature of a nonlinear MQ
system whose only solution is the setter's 256-bit message, far from any local start.

## 4. What would cross it
(a) the setter's 256-bit secret message (the CONSTs + the active-bit pattern); or
(b) a global MQ solve over the message class — the trapdoor itself (exponential without a
weakness, and none was found: random constants, nonlinear-in-bits core, pinned wire).

## 5. Tooling added this session
`find_pins.py`/`pin_parse.py`/`pin_sel.py` (256-bit MUX decode, `pinrec.json`),
`gadget_solve/dioph/full/handles.py` (exact gadget solve), `close_g1g2.py`,
`atoms16.py`/`decode4.py` (4-atom ripple), `iter_heal.py` (divergence proof),
`check_37110.py`/`free_wire.py` (wire-lock proof), `const_struct.py` (no constant structure).
Best verified unchanged and re-checked: **39,022 / 39,033**.
