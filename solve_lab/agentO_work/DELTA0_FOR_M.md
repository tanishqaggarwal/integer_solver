# delta0 — an exact lattice target for the 39,026 residual region

Agent O.  Everything is keyed by atom **source text** and **variable index** (shared
across models); atom *numbering* is not (E 40,727 / H 42,267 / M its own).  Equation
indices ARE shared and match `checker.py` line indices.

## What this is
The 39,026 witness residual lives in 9 atoms touched by exactly 13 equations, with 8
variables private to that region.  Over Q the system has a unique solution satisfying all
13; over Z five coordinates are blocked, by moduli `p, p, p, 2458959, 2458959*p`.
Solving `A z + B d = b0` over Z gives an integral boundary shift `d = delta0` supported
on exactly the four constants that are NOT p-multipliers.  Applying it makes all 13
region equations hold — **verified end-to-end in the region model**.

**It is NOT verified as realisable.**  Two of the four carriers are free variables and
cost nothing; the other two are derived, and their collateral is the open question.
That is what needs pricing.

## The four shifts

### atom `x_7068 - x_2099 - 7376877 * x_642`
- external expression: `x_7068 - x_2099`
- carrier variables: [7068, 2099]   (carrier factor 1)
- **shift = 2440 bits**, free to move by period exceeds 2458959*p (unmeasured)
- carrier increment (already divided by the carrier factor): 2440 bits
- external part of the atom that would define handle x_642; the private knob x_642 enters with coefficient -7376877, so this shift only matters modulo 7376877 = 3*2458959

### atom `x_4432 - x_19964 - x_28730`
- external expression: `x_4432 - x_19964`
- carrier variables: [4432, 19964]   (carrier factor 1)
- **shift = 2419 bits**, free to move by multiples of p
- carrier increment (already divided by the carrier factor): 2419 bits
- external part of the atom that would define handle x_28730; the private knob x_28730 enters with coefficient -1 here and +1 in the x_28730 - p*x_9413 atom, so only the DIFFERENCE of the two directions is a new degree of freedom

### atom `5113045 * (x_7075 * x_9118) - x_29854`
- external expression: `5113045 * x_9118`
- carrier variables: [9118]   (carrier factor 5113045)
- **shift = 2428 bits**, free to move by multiples of p
- carrier increment (already divided by the carrier factor): 2406 bits
- carried by the free variable x_9118 (zero-collateral in frame B); the shift may be moved by multiples of p, and gcd(5113045, p) = 1, so a representative divisible by 5113045 always exists

### atom `x_7075 * x_8731`
- external expression: `x_7075 * x_8731  (x_7075 = 1)`
- carrier variables: [8731]   (carrier factor 1)
- **shift = 2429 bits**, free to move by multiples of p
- carrier increment (already divided by the carrier factor): 2429 bits
- carried by the free variable x_8731 (zero-collateral in frame B)

## Which are free
- `x_8731` (atom `x_7075 * x_8731`) and `x_9118` (atom `5113045 * (x_7075 * x_9118) - x_29854`)
  are agent H's zero-collateral knobs in frame B — these two shifts cost nothing there.
- `x_7068 - x_2099` and `x_4432 - x_19964` are derived; **these two need pricing**.
  They are the external parts of the atoms that would define handles `x_642` and `x_28730`,
  i.e. two of the four handles the deliverable itself corrupts.

## Two simplifications worth using
- `x_642` enters `x_7068 - x_2099 - 7376877*x_642` with coefficient `-7376877`, and it is
  private, so the `x_7068 - x_2099` shift **only matters modulo 7376877 = 3 x 2458959**.
  That is a 23-bit condition, not a 2440-bit one.
- `x_28730` enters `x_4432 - x_19964 - x_28730` with coefficient `-1` and the
  `x_28730 - p*x_9413` atom with `+1`, both private, so only the DIFFERENCE of those two
  directions is a genuinely new degree of freedom.

## Do not scan configurations
Measured: the four boundary quantities are identically 0 across 35 configurations in E's
frame (1 distinct value out of 35 each) — a scan measures one point repeatedly — and the
admissible lattice has index >= 2^768 (hit rate ~2^-767).

Machine-readable copy: `DELTA0_FOR_M.json` (exact integers).
