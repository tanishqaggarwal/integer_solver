# RESUME_F — agent F (p-adic / multi-modular lifting)

## Best verified: 39,026 (shared baseline). My own best: **39,022** = `agentF_work/best_F_39022.json` (checker-verified).
## forward-closing the baseline's free inputs (`theirs_closed.json`). No improvement yet.

## Established (independently)
- checker baseline 39026 CONFIRMED (fails [12231,12270,12350,14584,18673,22044,29125]).
- Instance = straight-line circuit. 39,033 atoms; 30,001 definitions forming a **DAG**; 8,747 free inputs;
  **9,032 residual atoms**. All residual atoms = 0  =>  all equations hold.
- `fwd.py` Engine forward pass 16 ms; its equation score matches checker.py exactly.
- All-free-inputs-zero => **39,005 with only 3 broken residual atoms**.
- p = secp256k1 field prime, literal in the file (x26064). Handles enter as p*h. Prior claim CONFIRMED.
- The 3 atoms = "OR(a,b)=1" + "selected EC coordinate pair ≡ (K1,K2) mod p".
- Setting (a,b)=(1,1), coords (x22162,x30213)=(K1,K2), x24468=K1, x18956=K2 => **39,013**.
- Exact integer Jacobian: 5,747 affine knobs, 260 non-affine (254 OR-tree bits + 6 EC coordinates).
  Bipartite structure is near-diagonal; blocked components are size 8-9 and are inconsistent ONLY because
  the quadratic coordinate knobs are excluded from the linearization.

## Re-enter
    cd /home/user/integer_solver/solve_lab/agentF_work
    python3 parse3.py; python3 circ4.py; python3 sched.py; python3 supp.py   # rebuild (~1 min)
    python3 -c "import sys;sys.path.insert(0,'.');from fwd import Engine;E=Engine()"   # engine

## Next experiment (highest value)
Build the Jacobian at the BASELINE 39,026 free inputs (`theirs_closed.json`, 4 broken atoms
3130/3131/3132/7251), extract their bipartite components, and solve those components exactly over Z
INCLUDING the quadratic coordinate knobs (probe-based exact quadratic model, then Newton/HNF).
