# New instance (`EQUATIONS.txt` @ b616370) — progress & obstruction

Honest status: the re-randomized instance is the **same trapdoor family** but its twist is
**materially harder** than the previously-solved one. I reconstructed the partial
generically and isolated the twist, but have **not** cracked the final coupled system.

## What works (verified)
- Regenerated atoms (`poly_atoms.py`) and gates (`extract_atoms.py`): same family —
  degree histogram `{1:20124, 2:25480, 4:694}`, 841 huge-constant load atoms, a 220-var
  identity wire, ~694 perfect-square verifier atoms.
- **`rebuild_partial.py`** — a *generic* reconstruction (greedy topological gate
  orientation + forward-eval from free-inputs = 0, the same principle proven on the old
  instance where every free input is 0). Reaches **39,005 / 39,033** and isolates the twist
  to **4 nonzero atoms**. This confirms the core methodology transfers.

## The twist (4 atoms) and why it is harder
Two product-slack chains, exactly analogous to the old F/H:
- Chain 1: `x_14257 = x_7497·x_23917` must absorb `8863713·(x_18956 − BIGCONST)`.
- Chain 2: `x_32989 = x_11436·x_22399` must absorb the `HUGE2` gap.
- x_7497, x_11436 (and deeper x_22820, x_14393) are genuinely free (1 atom each).

**The decisive difference from the solved instance:** there the 220-var wire was a *free*
parameter, so I set it to 1 and the gaps absorbed as `1 · (−G)` — a trivial factorization.
Here the wire multipliers x_23917, x_22399 are **pinned to a fixed V0 ≈ 2²⁵⁶** by a hard
constraint (`x_26064 − V0` is a real shared atom, and a 5-hop identity chain ties
x_22399/x_23917 to x_26064). With the wire pinned:
- the product slacks can only produce **multiples of V0**;
- the gaps (BIGCONST, HUGE2) are **not** multiples of V0 (verified);
- so the non-V0 part of each gap must route into the two shared verifier checks
  **40907, 44255**, where the two chains' huge values must **cancel**.

The cancellation is not linear: `a·BIGCONST + b·HUGE2 (mod V0)` is ~2²⁵⁶ in both checks
(not reachable by the small linear vars), and the checks contain **products of
already-huge variables** (e.g. `x_29356·x_33469` with x_33469 ≈ 2²⁵⁶). Absorbing the
residue therefore requires solving a **nonlinear (product) Diophantine system** over ℤ —
the pinned wire has turned the old unit-factorization into a residue/factoring-type
problem. This is plausibly the intended hard core of this instance.

## Approaches tried (and why they stalled)
- `local_solve.py` — Smith-normal-form integer solve over spine+checks. Solves the spine,
  but the ℚ/linear treatment either blows up nonlinear check terms or ripples into
  neighboring atoms (39,005 → 38,963/38,992).
- `auto_solve.py` — auto-expanding frontier + integer solve. The coupled component grows
  without bound (the checks couple to thousands of vars), so a single SNF does not
  converge in memory/time.

## Honest assessment / next angles
The remaining core is a small but **nonlinear, V0-modular** integer system. Cracking it
likely needs one of: (a) the intended gate orientation that keeps the huge flow off the
small-coefficient checks (if one exists); (b) a targeted nonlinear solve that chooses the
product multipliers (x_29356, …) to hit the two check residues mod V0 while the V0-quotient
is handled by the free partners; or (c) recognition that this instance encodes a
factoring-hard verifier, in which case it is not solvable by local repair without the
trapdoor secret. Progress and tooling are committed; the 39,005/39,033 partial is real.


## DECISIVE UPDATE — the twist is global, not local

Setting the *natural* flow values (x_18956=BIGCONST, x_22820=BIGCONST//V0, the mod-V0
remainders, and the chain-2 analogs; see best/new_instance_partial_39007.json) zeros the
entire spine and reaches **39,007 / 39,033** — only the two verifier checks 40907/44255
(plus two resid-1 chain-3 atoms) remain. But the check residuals are ~2^300 and can only be
cancelled by products of already-huge variables (x_33469, x_24453 ~2^256), whose multipliers
then become huge and appear in *further* checks. Measuring the closure of this cancellation:

    closure from the 6 remaining atoms = 31,211 atoms / 26,915 variables

i.e. satisfying the checks requires re-deriving ~27k variables consistently — essentially
inverting the whole circuit. This is the intended trapdoor hardness. In the previously
solved instance the wire was FREE, which kept the twist LOCAL (224 vars); here the wire is
PINNED to V0, which makes the same twist GLOBAL. Local repair / linear solving cannot reach
it, and the exact linear systems blow up (entries ~2^256 compound through elimination).
Conclusion: this instance is solvable only with the setter's message (the consistent global
assignment), not by the local-construction method that cracked the free-wire instance.

## FRESH FIRST-PRINCIPLES ANALYSIS (no prior framing)

Treating the system as a raw polynomial object:
- Linear subsystem rank 19,406 / 28,837 vars; iterated linear+constant closure pins only
  6,683 / 38,748 vars — NO linear collapse.
- 8,581 "input" (never-gate-output) vars; they carry **0 independent linear constraints**
  (Schur complement of the linear system onto the inputs is empty).
- Every atom LINEARIZES by substituting gate outputs for products, EXCEPT the 768 perfect-
  square verifier atoms. So the entire nonlinearity = product-consistency (auto-satisfied by
  forward-eval) + the 768 verifier squares.
- => the problem is exactly: choose the 8,581 free inputs so the 768 verifier squares hold.
  This is circuit inversion; there is no linear leverage on the inputs.
- Irreducible obstruction (independently re-derived): at inputs=0 all verifier squares hold
  and only the LINEAR gap atoms fail; satisfying a gap forces an input to carry
  BIGCONST mod V0 ≈ 2^256 (since V0 ∤ BIGCONST), and that same input feeds the verifier,
  breaking it. Linear gaps and nonlinear verifier pull the inputs oppositely; the V0-remainder
  is unabsorbable. Confirms the pinned-wire hardness from pure algebra.

## BREAKTHROUGH (Session, 7hr mandate): twist mechanism cracked via circuit inversion

The pessimistic "global/unsolvable" conclusion was WRONG. Working the raw equations:
- The 2 verifier checks (40907/44255) are each a linear combo of the SAME 17-gate sub-
  circuit; they vanish iff those gates hold. The gaps route into FACTORABLE product gates
  (x_38045=x_15298*x_22162, x_10156=x_15298*x_30213).
- Recursive circuit inversion (achieve2.py): activate x_15298=1 (boolean OR-gate tree), set
  value inputs x_30213=BIGCONST, x_22162=H2. This makes x_37892=BIGCONST, x_13682=H2 and
  FIXES THE ENTIRE MAIN TWIST — gaps 602/1465 AND both hard checks 40907/44255. Verified:
  those atoms go to 0.

## Remaining core: the 256-bit message codeword
- x_15298's activation tree has 256 free-input leaves, ALL 'load bits' (free inputs with a
  >10^40 coefficient). Activating x_15298 REQUIRES setting >=2 load bits (one per OR-tree
  for x_7715, x_34554); each triggers a huge load (HUGE*bit).
- Each load atom (coef*x_A + HUGE*bit - bit*absorber) is absorbed by a FREE-INPUT absorber
  (x_38460, x_27466, ...). BUT those absorbers reappear in verifier squares (e.g. atom
  18081 = x_38460*x_c, needing x_c=0), coupling back into the deg-2 verifier.
- => the 256 load bits + their absorbers must jointly satisfy the verifier squares. This is
  the genuine trapdoor core ("256-bit core"): the bits load huge values (linear in bits),
  the verifier checks them (quadratic). Best full partial remains 39007/39033 (natural.json);
  the inversion demonstrably solves the twist but the codeword search remains.

## Creative reduction: single-clean-bit routing (12 -> 5 broken)
- Insight: x_15298=x_7715*x_34554 is NOT needed. Since x_23597=1 at base, x_34606=x_7715, so
  BOTH gaps route through x_7715's products alone: x_37892=x_7715*x_16742 (x_16742=BIGCONST),
  x_13682=x_7715*x_12186 (x_12186=H2). Activating just x_7715 halves the load-bit cascade.
- Forcing activation through a CLEAN bit (x_36314, whose absorbers avoid the deg-4 verifier
  squares) + achieve3/4/5.py reduces the broken set to 5 atoms (39004/39033).
- IRREDUCIBLE COUPLING (the MQ): the clean bit's load absorbers (e.g. x_14282, x_34734) are
  SHARED between their load atom AND a verifier check (atom 42167), which also carries the
  twist's huge gap-vars. Each absorber is over-determined: it must both absorb the bit's load
  AND satisfy the verifier. Consistency requires the RIGHT bit pattern = the 256-bit codeword.
- Even the "cleanest" bits have absorbers in 7+ atoms; no isolated activation exists. So the
  reduction bottoms out at the same MQ core, now sharply localized: choose activation bits so
  that {load-absorption} and {verifier checks} agree. Best partial remains 39007/39033.

## DEFINITIVE STRUCTURAL MAP (this session) — the core is a 256-bit load-absorption message

Re-derived the whole instance from the raw equations, treating it as a fresh object:

1. **The gate DAG is fully ACYCLIC** (Tarjan SCC: 36,418 vars, 0 cycles, 0 self-loops).
   Every gadget has canonical form `(output − f(inputs))`, so orientation is *readable*, not
   guessed. Forward-eval from any free-input choice satisfies **every wiring equation
   automatically** (each gadget = 0 ⇒ each E = 0 ⇒ c·E=0 and E²=0 hold).

2. **The only obstruction is the 1,841 join points** (vars defined by ≥2 gadgets that must
   agree) plus square-only vars. Localizing the 28 all-zero failures: culprits are
   `x_2300` (=1 pinned AND =x_9274 ⇒ **x_9274=OR(x_7715,x_34554)=1 forced** — why all-zero
   fails), `x_24468=C1` & `x_18956=C2` (huge ~2^296 forced constants, C1/C2 in huge_consts.json).

3. **The twist is a MUX** on two control bits (x_7715,x_34554). Quadrant selectors
   x_15298=AND, x_34606=x_7715·(1−x_34554), x_5647=x_34554·(1−x_7715). Product gates route
   data. With control **(1,0)**: x_37892 = x_16742, x_13682 = x_12186 (free inputs route
   DIRECTLY to the two forced constants — clean, no huge products). So set x_16742=C2,
   x_12186=C1. Symmetric routes exist for (0,1) via x_24908,x_14853 and (1,1) via x_30213,x_22162.

4. **The 256-bit message.** The spine's backward cone (2,411 vars) contains exactly **256 free
   boolean bits** = the *disjoint* union of x_7715's cone (178 bits) and x_34554's cone (78
   bits); the two cones share 0 inputs. All 178 bits single-flip-activate x_7715 (it is
   effectively a giant OR). The bits reach the spine ONLY through the 2 control values.

5. **Why activation isn't free — LOADS.** Each message bit appears in verifier squares of the
   form (e.g. eq84) `x_a·x_b = bit·(x_load − HUGE)`. With bit=0 the square is satisfied by
   zeroing a product (trivial — this is why all-zero passes all squares). With bit=1 it forces
   a **huge load** x_load ≈ HUGE that must be absorbed downstream. Every quadrant's data input
   is shared across 31–48 squares, so setting it perturbs them all.

6. **Measured hardness.** Best single-bit activation = 39,007/39,033 (x_22106), same as the
   natural partial. Fixing the activation bits and solving the residual failing equations as a
   **linear system in the free data/absorber inputs is INCONSISTENT** (SNF: no integer
   solution) — proving that no data choice repairs a wrong message; the 256 bits must be
   exactly the setter's codeword. Confirmed genuine MQ/codeword trapdoor (matches the
   pinned-wire hardness: here the twist multiplier wire is pinned to V0, unlike the solved
   free-wire instance).

Tools added: scc.py, localize.py, coupling.py, build_twist.py, scan_bits.py, newton.py,
examine_fail.py. Best verified partial remains best/new_instance_partial_39007.json (39,007).
Next probes: smaller x_34554 activation region (78 bits, quadrant (0,1)); multi-bit
load-cancelling combinations; per-bit load-target absorbability ranking.
