# Agent D log

- t0: verified `best/new_instance_partial_39026.json` = 39,026 / 39,033, failing
  [12231, 12270, 12350, 14584, 18673, 22044, 29125]. Prior claim CONFIRMED.
- built own caches (build_cache.py): 42,267 atoms, 39,033 eqs, 31,475 gates,
  10,792 checks, 7,273 free inputs, topo covers 29,675/31,475 (1,800 in cycles).
- forward eval from witness free inputs -> 38,996 (6 nonzero checks). Confirmed prior.
- zero x_9118,x_8731,x_1329,x_10903 -> 39,004. + force x_24548:=x_25442,
  x_14853:=x_1308 -> 39,002 fixed point, residual = {a19297,a19299,a21617,a30984,
  a36185,a37662,a40812}. (D_state1.json)
- engine2.St: incremental apply/revert, ~5 ms/move exact.
- rad.py reverse AD mod 2^61-1: knob list per atom in 0.02 s.
- big-constant census: only 4 public constants in the instance (p, C_A=x_24453,
  C_B via a688, C_C via a1618).
- decompiled a688/a1618: they pin y3 and x3 of the EC addition to C_B/C_C mod p.

## Breakthrough chain (session part 2)
- adv3.py (advice Gauss-Seidel over the 192 congruences u ≡ w mod p) : D_state1 -> **39,013**
  (D_adv.json). Residual = {a19297,a19299,a30984} + 2 one-equation shadows, and these are
  exactly `sel * (combination of A,B) ≡ 0 (mod p)` with A = x_35389, B = x_6671.
- Decompiled exactly: A = (x1+x2+x3+C_A)(x2-x1)^2 - (y2-y1)^2 ; B = (y3+y1)(x2-x1) - (y2-y1)(x1-x3);
  verified digit-for-digit against x_35389 / x_6671, and x_9192 = x1+x2+x3 verified numerically.
- scanAB.py: all 7,273 free inputs perturbed + advice re-solve. ZERO cost-free movers of (A,B).
  Cheapest: x_22162 (=x3) cost 2, x_30213 (=y3) cost 4, then x_22152/x_33462/x_6418/x_12553 at 10-14.
- ecsolve2.py: A,B are EXACTLY affine over Z in (dx3, dy3); solved
  p | x_11150, 6672769*p | x_25739, p | x_37758 by CRT -> then closed the three handles
  -> **39,017** (D_39017.json), only 3 nonzero atoms {a688, a1618, a40608} = the x3/y3 pins.
- table.py/condpins.py: 256 free boolean bits, each pinning EXACTLY 2 vars to 2 296-bit
  constants; 178 bits gate x_7715 (P1), 78 gate x_34554 (P2). bitswap.py CONFIRMS
  (x1,y1) = the selected bank-1 pair, (x2,y2) = the selected bank-2 pair.
- banks.py: exhaustive 178 x 78 x 4 orderings -> NO pair gives A=0 or B=0. The one-bit-per-bank
  branch is unsatisfiable; the constants are NOT secp256k1 points (0/256 on curve, 278/512 QR).
- two bits on in one bank -> x1 = 0 (not a subset sum; it is a tree MUX).
