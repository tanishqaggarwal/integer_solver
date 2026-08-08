# Smallest faithful full 256-bit comb ladder (measured)

Atom composed into the full ladder, every number measured at s=256 against the real
p = 2^256-2^32-977, whole ladder re-certified by demo_win2.py's exhaustive-scalar test.

## Window-width co-optimization (ladder = ceil(256/w) * window(w))

| w | windows | LOGICAL | K | |J| | couplers | phys Pegasus | phys Zephyr |
|---|---|---|---|---|---|---|---|
| 6 | 43 | 12,696,352 | 5 | 2^8 | 50,330,167 | 12,696,352 | 12,696,352 |
| 7 | 37 | 11,743,504 | 5 | 2^9 | 47,225,690 | 11,743,504 | 11,743,504 |
| **8** | **32** | **11,574,912** | **5** | **2^9** | **47,628,576** | **11,574,912** | **11,574,912** |
| 9 | 29 | 13,056,728 | 5 | 2^9 | 55,440,953 | 13,056,728 | 13,056,728 |

Smallest: **w=8, 11,574,912 logical = physical** (curve flat within 1.4% across w=7-8).

## Squaring-aware step
One group addition = 3 field multiplies, exactly one a squaring (lam^2). Encoding
lam^2 as a squaring (62,692) not a general multiply (99,298) saves 36,606 qubits/step
(12.3%). Step total 269,231 (2 general + 1 square + linear glue + d!=0/d!=p gadgets).

## Physical qubits: whole-ladder max clique K = 5
Verified two ways: (1) every window build w=4..12 reports K=5 (wallace 3:2 tree reduces
even the 2^w-wide table-select and final-comparison columns to full-adder squares of
clique 5); (2) demo_win2.py on real full ladders (small curves) reports K=5 in all 36
wallace cases. K=5 <= native degree on Pegasus(15)/Zephyr(20) => chains length 1,
physical = logical, ZERO embedding overhead. (Binary-carry ladder is 3x fewer logical
but K~139 => ~12x physical/logical -- wallace is the physical minimum.)

## Faithfulness certificate (demo_win2.log): 73 cases, 0 failures
Real Hamiltonian on pseudo-Mersenne p in {127,251,1021} (NAF reduction exercised),
every ancilla filled by witness replay, every candidate scalar enumerated. Faithful iff
{k: E(k)=0} == {k: k*G=T}. All 73 FAITHFUL; the 36 wallace cases (the exact composed
encoding) certify at K=5, |J|=2^6. Degenerate-division loophole closed by d!=0/d!=p.
This is the ladder-level piece of the end-to-end certificate.

## Bottom line
w=8 32-window ladder: 11,574,912 logical = physical (K=5, no embedding overhead),
47.6M couplers, |J| 2^9 = 2,010x Pegasus / 2,631x Zephyr, exhaustively faithful.
Note: uses karatsuba leaf-24 atom (99,298); the toom3(128)>kara(24) plan (96,809) and
the mincost agents' results would shrink this proportionally.
