# Agent G — RESUME

## Best verified score
**39,026 / 39,033** — the inherited `solve_lab/best/new_instance_partial_39026.json`,
re-verified: `python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
-> `satisfied 39026/39033 (7 failing)`, failing `[12231,12270,12350,14584,18673,22044,29125]`.
I did not beat it; nothing was written to `best/`.

## WHAT THE INSTANCE IS (settled, exactly, this session)
The circuit is a **secp256k1 double-and-add**. `secp_identification.json` holds every
constant; `g41_curve.py`, `g50_points.py`, `g51_chain.py`, `g53_export.py` reproduce it.
* Curve in the instance: `y^2 = x^3 + a2 x^2 + a4 x + a6` with
  `a2 = K = 97553848499418123410591666447050222001188385549510401465815187079080512838891`
  (K is exactly the constant in my A identity), `a4 = a2^2/3` so **j = 0**,
  `a6 = 77755683306591771556999954628254672912734268662742093169295805431582354953490`.
  After `x -> x - a2/3`: A_short = 0, `B_short = 640195336800308764084431987622108290587517006345542821859873258203935985247 94`,
  and `B_short/7` is a 6th power -> **isomorphic to secp256k1 over F_p**.
  Map: `x_sec=(x+a2/3)/u^2`, `y_sec=y/u^3`,
  `u = 4210889811980686189396764679825672592540066047176031544704936155054310740018`.
* **Exactly 256 boolean free inputs carry huge load pins, and their points form a perfect
  doubling chain: P(bit at chain index i) = [2^i]P0 for all 256 positions, 0 exceptions.**
  P0 has order n. `P1 = [2^72]P0` (bit x24601), `P2 = [2^235]P0` (bit x2081) — the two bits
  set in the base state. P3 is the pinned target, on the curve, order n.
* My exact reduction shows the whole instance mod p is `A = 0 and B = 0`, i.e.
  **P3 = P1 + P2**; with m message bits on this becomes the m-term double-and-add.
* **A full solve is exactly `k = log_{P0}(P3)` on secp256k1.** The bits act by group
  doubling, NOT affinely on (A,B) — so no LLL / low-density subset-sum attack applies
  (`sum b_i [2^i]P0 = P3` *is* the discrete log; density 1, group exponential map).
* Hamming weight of k <= 4 is ruled out by meet-in-the-middle (`g52_lowweight.py`).

## The reusable machinery (exact, validated)
`gsym.py` / `gsym2.py`: exact symbolic forward evaluation over F_p (legal because every
gate output coefficient is +-1). With ALL 6,117 non-boolean free inputs symbolic the pass
takes 0.5 s, 0 gates skipped, and of the 10,792 check atoms only **2,029 are non-constant**
(1,883 linear / 141 quadratic / 5 cubic). Sparse F_p elimination (`gsolve.py`): rank 1470,
4,647 free params, consistent; substitution turns every nonlinear check into a CONSTANT,
five nonzero. Validated at random points: 0 mismatches on all 10,792 checks; and at the
39,026 deliverable's own point, 0 atom mismatches over all 42,267 atoms with exactly 7
nonzero equations. With the 112-symbol closure the system is 57 polynomials / 196 monomials.
Equation-level version: 6,774 non-trivial equations, forcing all 6,613 linear ones leaves
exactly 20 nonzero = AG_39013's 20 failing.
**Every rank/kernel/ceiling in this repository before this was a tangent-space quantity;
these are the exact polynomial system.**

## Re-enter
```
cd /home/user/integer_solver/solve_lab/agentG_work
python3 g41_curve.py           # curve identification from my own A,B identities
python3 g51_chain.py           # the 256-point doubling chain and its root P0
python3 g53_export.py          # -> secp_identification.json (re-verifies 256/256)
python3 g23_allsym.py ; python3 g24_bigsolve.py base ; python3 g35_eqsolve.py -
python3 g29_frame.py - 2081 24601 4287 13195      # exact reduce per boolean frame
```
NOTE `s9/eff/lib.py` does `os.chdir(solve_lab/s9)`; use absolute output paths.

## Highest-value next experiment
There is no cheap algebraic route left: the residual is the discrete log. The only honest
options are (a) extend `g52_lowweight.py` to Hamming weight 8-10 by meet-in-the-middle
over the two 128-bit halves (~1e7 hash entries, hours, cheap insurance in case the setter
picked a sparse k), or (b) accept 39,026 and spend the remaining effort on the
minimum-weight coset-decoding problem in my exact equation-level model (6,613 linear +
161 nonlinear equations, 4,652 unknowns) to see whether 7 failing equations can be beaten.
