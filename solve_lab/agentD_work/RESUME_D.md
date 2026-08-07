# Agent D — resume brief

## Best verified so far
Global best is still `solve_lab/best/new_instance_partial_39026.json` = **39,026** (re-verified).
MY best own state: `agentD_work/D_39017.json` = **39,017**, checker-verified,
failing [56,133,2071,8073,13660,15299,16622,17726,19066,22093,25480,28653,31061,32894,34517,34892].
Its residual is only THREE atoms {a688, a1618, a40608} = the two mod-p pins on (x3,y3).
Pipeline that produced it: D_state1 -> adv3.sweep -> 39,013 (D_adv.json) -> ecsolve2.py -> handles -> 39,017.

## Pipeline (all inside agentD_work/, independent of s9/s10 caches)
```
cd solve_lab/agentD_work
python3 build_cache.py          # rebuilds cache/{atoms,polys,gates,topo}.pkl  (~45 s)
python3 -c "import dlib as L, engine2 as E; st=E.St(L.load('D_state1.json')); print(st.score, st.nz())"
```
* `dlib.py`   — atoms/polys/gates/topo + exact eval (42,267 atoms, 39,033 eqs, 31,475 gates, 7,273 free inputs)
* `engine2.py`— `St` state with **incremental apply/revert** (~5 ms per move, exact ints)
* `rad.py`    — reverse-mode AD mod 2^61-1 → free-input knob list for any atom in ~0.02 s
* `fwd.py`, `dig.py`, `ortree.py`, `consts.py`, `pins.py` — structure tools

## Established by me this session (measured, not inherited)
1. Forward eval from the 39,026 witness's free inputs → 38,996, 6 nonzero checks. Confirmed.
2. Zeroing free inputs x_9118,x_8731,x_1329,x_10903 → **39,004**; then forcing
   x_24548:=x_25442, x_14853:=x_1308 → fixed point **39,002** with residual moved to
   {a19297,a19299,a21617,a30984,a36185,a37662,a40812}. (Prior "oscillation" confirmed.)
3. **Only 4 large public constants in the whole instance**: x_26064 = p (secp256k1),
   x_24453 = C_A (256-bit), a688: x_18956 ≡ C_B (296-bit) mod p, a1618: x_24468 ≡ C_C (296-bit) mod p.
4. Decompiled the pins: with selectors x_15298=1, x_34606=x_5647=0,
   **x_18956 ≡ y3 and x_24468 ≡ x3**, i.e. a688/a1618 pin the *output point* (x3,y3) = (C_C, C_B) mod p.
   In D_state1 x_22162 = C_C and x_30213 = C_B exactly, so those pins already hold.
5. x_15298 = AND of two 192-leaf OR trees, 384 leaves, 256 of them free; only x_2081 and
   x_24601 are nonzero. Zeroing both → 38,879 (selector door is expensive). Confirmed prior.
6. Residual conditions at D_state1 reduce to **A ≡ B ≡ 0 (mod p)** (EC addition identities)
   plus the advice equality a21617: x_14623 ≡ x_27522 (mod p).

## Next experiment (highest value)
Solve A=B=0 in closed form: pick x1,x2 free, s = sqrt(x1+x2+x3) mod p,
y1 = s(x1-x3) - y3, y2 = y1 + s(x2-x1) — then re-solve the advice equality DAG and measure.
The advice checks are `K*(u - w) - p*handle` with u a free input, so each is closable by
u := w; the question is whether the DAG fixed point is compatible with the EC solution.

## Structure, fully decompiled (mine, measured)
* 256 free boolean bits; each pins EXACTLY 2 variables to 2 296-bit constants (condpins.py).
  178 gate the OR-tree x_7715 ("P1 selected"), 78 gate x_34554 ("P2 selected").
  bitswap.py: turning bit i on (and setting its 2 pinned vars) makes (x1,y1) = that bit's pair.
* Selectors are mutually exclusive: x_15298 = f1*f2, x_34606 = f1*(1-f2), x_5647 = f2*(1-f1),
  x_15574 = 1-x_15298.  With x_15298=1 the two public pins read  x3 ≡ C_C,  y3 ≡ C_B (mod p).
* Whole instance is satisfied iff (at the advice fixed point) A = x_35389 ≡ 0 and B = x_6671 ≡ 0 (mod p),
  A = (x1+x2+x3+C_A)(x2-x1)^2 - (y2-y1)^2,  B = (y3+y1)(x2-x1) - (y2-y1)(x1-x3).
* A,B are EXACTLY affine over Z in (δx3, δy3); solving them costs breaking the two pins = 16 equations.
* banks.py: exhaustive over 178x78x4 orderings — NO bit pair makes A or B vanish. Table constants
  are NOT secp256k1 points. Two bits on in a bank gives x1 = 0 (tree MUX, not a subset sum).

## Next experiment
Rank ALL knob pairs (u1,u2) that span (dA,dB) mod p by the number of equations their broken
atoms occupy, and run the ecsolve2 CRT construction with the cheapest pair. Currently
(x3,y3) costs 16; a pair costing <= 6 would beat 39,026.
