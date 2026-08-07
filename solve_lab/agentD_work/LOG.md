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
