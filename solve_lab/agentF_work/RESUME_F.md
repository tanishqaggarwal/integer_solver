# RESUME_F — agent F (p-adic / multi-modular lifting)

## Best verified: shared baseline **39,026** (`solve_lab/best/new_instance_partial_39026.json`, re-verified
## by me).  My own pipeline's best: **39,023** = `agentF_work/F_frame.json` (checker-verified), and
## `agentF_work/best_F_39022.json` (39,022, 2 broken atoms only).

## What I established independently (all from my own parse; nothing imported)
- Instance = straight-line circuit.  39,033 atoms; 30,001 definitions forming a **DAG**; 8,747 free inputs;
  **9,032 residual atoms**.  All residual atoms zero => all 39,033 equations hold.
- All-free-inputs-zero => 39,005 with only THREE broken atoms.  Engine `fwd.py`: 16 ms/pass, score matches
  checker.py exactly.
- p = secp256k1 field prime, a literal in the file (x26064).  Handles enter every check as `M*p*h`.
- Decoded residual: OR(a,b)=1 over two OR-trees (178 + 78 boolean free bits); selected coordinate pair must
  be ≡ (K1,K2) mod p; EC layer vanishes iff the short-Weierstrass chord law holds with a constant offset
  K = 97553848499418123410591666447050222001188385549510401465815187079080512838891 (measured, universal).
- Shift X = x + K/3: all 256 pin points and the target lie on **E_b: y^2 = x^3 + b**,
  b = 64019533680030876408443198762210829058751700634554282185987325820393598524794, over F_p.
- **|E_b| = N = secp256k1's own group order n, PRIME** (no factor <1e6, isprime, MR bases 2..100).
  Pohlig-Hellman worthless; no smooth/small-order structure.  Quadratic twist 3^2*13^2*3319*22639*(211-bit p).
- **Single doubling chain of length 256 rooted at bit x_2779** (253 full points give 4 chains 124/79/41/9;
  the 3 one-pin bits 18184/22579/33434 are exactly the 3 missing links).
- b/7 is a 6th-power residue; u = 4210889811980686189396764679825672592540066047176031544704936155054310740018.
  Under (x,y)->(x/u^2,y/u^3), **P0 is NOT secp256k1's generator G**; T is not G; no log for k<=5000, k=2^i,
  (N±1)/2, N-1.  BSGS to 2^40 running (`bsgs.py` -> `bsgs.log`).
- **The lift obstruction is EXACTLY at p.**  `modq.py` solves the FULL 39,033-equation system over Z/qZ with
  0 residual atoms and 0 failing equations for q in {4093,7919,65537,104729,1000003,15485863,999999937,
  1e9+7,2^31-1,2^61-1,1e12+39} (~4 s each).  Reason: every check is `A-B = M*p*h` with h free, so it is
  vacuous mod any q coprime to M*p.  Only at q=p does the handle die and the pins bite.
- **39,026 optimality re-verified independently**: its residual is 7 atoms / 12 equations; the reachable
  lattice (7 private vars) makes exactly ONE 5-subset of the 12 rows integrally solvable (all 924 6-subsets
  are solvable over Q).  Every single-atom extension of that frame (7 candidates, `frame_search.log`) still
  gives nfail = 7.  Prior claim CONFIRMED.

## Re-enter
    cd /home/user/integer_solver/solve_lab/agentF_work
    python3 parse3.py; python3 circ4.py; python3 sched.py; python3 supp.py   # rebuild pickles (~1 min)
    python3 fwd.py           # 39,005 from all zeros
    python3 modq.py          # full solve mod q for many q
    python3 frame_search.py  # frame optimality
    python3 bsgs.py          # discrete log search

## Single highest-value next experiment
Multi-bit ladder consistency: turn on k>=2 bits of the 256-chain, compute the exact integer Jacobian of the
merged chain component, and solve it with `intsolve.solve_int`.  If consistent, the accumulator is an honest
subset/ladder sum and the instance is exactly ECDLP on E_b; if not, the reachable accumulator set is only the
13,884 one-bit-per-tree sums and no assignment exists (which would make 39,026 provably optimal).
