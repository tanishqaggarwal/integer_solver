# RESUME_E — agent E private checkpoint  (keep <=30 lines)

**WITHDRAWN ANGLE:** generator/authorship reverse-engineering (PRNG hunting, emission-order
forensics, "run the generator forward").  Do NOT resume it.  Everything below is computed from
the equations as mathematical objects (syntactic sub-expression decomposition + propagation).

## Established (all reproducible, exact integer arithmetic)
1. `best/new_instance_partial_39026.json` verifies at 39026/39033.  CONFIRMED.
2. `parse3.py` -> `model3.pkl`: every equation = outer_scalar * Z-linear combo of **atoms**
   (9,710 equations are perfect squares S*S).  40,727 distinct atoms / 38,748 vars.
3. `aeval.py <assign>`: at the 39,026 deliverable exactly **8 atoms** are nonzero and they
   reproduce the failing set exactly -> the atom model is faithful.
4. `dag.py` -> `dag.pkl`: 35,004 atoms are definitions `x_out - RHS`; def relation is ACYCLIC;
   8,365 vars never defined ("free").
5. **`harness.py`** = fast forward map free-inputs -> all 38,748 vars -> violated atoms
   (~0.14 s per full evaluation; orientation cached in `orient.pkl`).
6. **free=0 leaves only 3 violated atoms** -> 38,998/39,033 (`prop2.json`).
   Seeding `x_18956 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626`
   -> violated atoms {20212, 20215, 24403} -> **39,009/39,033** (`prop_s1.json`).
7. **Decoded residual** (this is the whole remaining obstruction in this frame):
   - a24403: `1 = x_29237 - x_23134` with `x_29237=a+b`, `x_23134=a*b`,
     `a=x_7715=OR(x_8599,x_21839)`, `b=x_34554=OR(x_25956,x_7304)`  => **OR(a,b)=1 forced**
     (x_9274=x_2300=1 is pinned by a literal).  Currently a=b=0.
   - a20212: `x_13913 = a*(1-b)*x_12186 + b*(1-a)*x_14853`  (a 2-way MUX; x_12186, x_14853 free)
     with `x_13913 = x_13682 - x_15298*x_22162`, `x_15298 = a*b`.
   - a20215: `x_24530 = x_5647*x_24908`, `x_5647 = b*(1-a)`.

## Best verified so far: 39,026 (the pre-existing deliverable).  My own best: 39,009.
## Next experiment
Force `a=1` (or `b=1`) by satisfying the isZero gadget upstream of x_8599 / x_21839 /
x_25956 / x_7304, then re-propagate; then set free x_12186 (resp. x_14853) equal to the
required MUX output.  Use `harness.forward(seed_dict)` + `harness.eqfails`.
