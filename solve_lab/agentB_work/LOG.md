# Agent B — RESUME (independent re-parse track)

## Established so far
- Verified the lab claim myself: `best/new_instance_partial_39026.json` -> 39026/39033,
  failing [12231,12270,12350,14584,18673,22044,29125]. CONFIRMED.
- Built an INDEPENDENT parser + model of EQUATIONS.txt (no lab code inherited):
  - `bparse.py`  : recursive-descent parser (v1, flattens '-' into *-1)
  - `bmodel.py`  : expanded-atom model     -> model.pkl  (94,202 atoms)
  - `bmodel2.py` : factored-atom model     -> model2.pkl (91,365 atoms, 87,083 factors)
  - `bdag.py`    : hash-consed DAG         -> dag.pkl (1.16M nodes)
  - `bparse2.py`/`bmodel4.py`: gate-faithful v4 -> model4.pkl (40,721 atoms, 42,553 factors)
  - `bparse3.py`: v5 parser that KEEPS paren groups ('g' nodes)  <-- current work
- STRUCTURE (my own derivation, matches nothing inherited):
  * every equation is `scalar * L` with L a linear form over gate atoms;
    kinds: sL 18478, pow2 8927, same 6589, plain 4256, pow4 783. No exceptions.
  * modulus p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
    (= secp256k1 p = 2^256-2^32-977) appears as a literal; 516 literals of 287-296 bits.
  * KEY PARSE FINDING: the generator wraps the HEAD gate of each packing chain in an
    EXTRA paren group -> that is the reliable gate boundary signal (used by bparse3.py).
- Naive "zero every gate" propagation (bprop.py) scores 26,915/39,033. Not competitive yet.

## Best artifact so far
- Nothing above 39,026 yet. Reference: /home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json

## Re-enter
    cd /home/user/integer_solver/solve_lab/agentB_work
    python3 bmodel4.py           # ~30 s, rebuilds model4.pkl from scratch
    python3 bgates.py            # gate/definition census

## Next experiment
Finish `bmodel5.py` (uses bparse3.flat_pack) -> exact gate decomposition, then a
topological forward evaluation of the circuit from free inputs, checking whether
ALL gates can be zeroed simultaneously.
