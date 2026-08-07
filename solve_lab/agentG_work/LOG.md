# Agent G LOG

## t0 — orientation
- Read PROMPT.txt, RESUME.md (1129 lines), STATE.json.
- Verified `best/new_instance_partial_39026.json` -> 39026/39033, 7 failing. CONFIRMED.
- Env: python-flint 0.9.0, sympy, z3, pysat, numpy. No scipy/sage/msolve/networkx/gmpy2.

## Exp G01 — coordinate/curve probe (`g01_probe.py`)
Prior claim (Part XXVI) says the core is secp256k1 point addition with
x1=x12186,y1=x16742,x2=x14853,y2=x24908,x3=x22162,y3=x30213.
Measured at 4 states:
- 39026 deliverable: x1==x2, y1==y2 exactly (mod p) -> the "addition" is degenerate there.
- AG_39013 / EC_39014 / PF_39015: P1 != P2.
- **None of P1,P2,P3 satisfies y^2 = x^3 + 7 (mod p) at any state.**
=> The EC reading is a formal algebraic identity, not literal secp256k1 points.
   Do NOT expect curve-theoretic structure (order n, discrete logs) to apply.
