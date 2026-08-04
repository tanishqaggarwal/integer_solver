# SOLVED — complete integer solution to EQUATIONS.txt

**All 39,031 / 39,031 equations satisfied exactly in ℤ.**
Solution file: `best/SOLUTION.json` (also at repo root `SOLUTION.json`), 38,748 vars, 2,954 nonzero.
Verify: `python3 checker.py best/SOLUTION.json` → `RESULT: OK`.
Independently re-verified by re-parsing `EQUATIONS.txt` (distinct code path): 39031/39031, 0 failing.

## The winning move (Session 8)

The previous sessions worked in the **confluent forward-eval orientation** and the **mod-P
reformulation**, where the twist looked like a rigid collision `x_3183=x_17728 & x_9770=x_18274`
requiring control bits (`x_12779≥2`), a div-wire escape, and a high-dimensional co-activation of
thousands of free vars. All of that was an artifact of the frame. Working directly in the **raw
equation space** collapses it.

### 1. At best_partial (39,019/39,031) only TWO atoms are nonzero
Evaluating the gate residuals at `best_partial_39019.json`, every atom is 0 except:
- **H** = `(x_17728 − x_3183) + x_9982`  (with x_9982 = 0, so H = G, the invariant gap)
- **F** = `6033033·(x_18274 − x_9770) − x_26977`  (with x_26977 = 0)

All 12 failing equations are linear combinations of the shared atoms that include H or F.

### 2. The gap is NOT rigid — it is a product slack
In the RAW equations, `(x_17728)-(x_3183)` **never** appears alone: all 16 occurrences are
`((x_17728)-(x_3183)) + x_9982`. So the "rigid a44271: x_3183=x_17728" was a reformulation
artifact. The real constraints are two **product-slack** activations:
- H = 0  ⟺  `x_9982 = −(x_17728−x_3183) = −G`, and atom 1818 sets `x_9982 = x_12518·x_9897`.
- F = 0  ⟺  `x_26977 = 6033033·(x_18274−x_9770) = F0`, and atom 1816 sets `x_26977 = x_20510·x_31302`.

`G  = 63398753350954830538284979531311478224817569395477016427713014637060524103217265241016814`
`F0 = −167517014178196647187119670838426778498383583571049457784393352496726267704305162856182499354817`

### 3. The hubs live in a QUIET 220-var wire
`x_12518` and `x_20510` are hub variables (271 / 237 equations) — but both belong to the **same
220-variable identity class** as `x_15` (the "wire"; members satisfy `x = ±x_15` via 2-term
identity atoms). Setting the whole wire uniformly to `sign·V`, holding all non-wire vars at
best_partial, leaves **every** atom's V¹ and V²⁺ coefficient **exactly zero** (integer, all 5233
touched atoms). The wire is a genuinely free parameter, decoupled from best_partial.

### 4. Direct construction (V = 1)
- Move the wire: every member → its sign (`x_12518 = x_20510 = +1`).
- Partners (each appears in only 2 atoms — its slack-def and one verifier square):
  `x_9897 = −G`,  `x_31302 = F0`.
- Slack outputs: `x_9982 = −G`,  `x_26977 = F0`.
- Everything else stays at best_partial.

This satisfies atoms 1818, 1816, H, F — and, as the exact checker confirms, **all 39,031
equations**. The verifier square atom 40782 = Q² (Q a 38-term deg-2 form) also lands at Q = 0.
`x_12779` and `x_24026` stay 0 throughout — the div-wire / dirty-bits path was never needed.

## Reproduce
`solve_lab/build_solution.py` regenerates `best/SOLUTION.json` from `best_partial_39019.json`
in one pass and runs the checker.
