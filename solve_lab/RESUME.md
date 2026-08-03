# RESUME — read me first

## Status
Best verified partial: **39,019 / 39,031** equations (exact in ℤ).
File: `best/best_partial_39019.json`. Verify: `python3 checker.py best/best_partial_39019.json`.

## THE problem is now fully reduced (Session 6)
The entire remaining obstruction (4 atoms: 1817, 30378, 40782, 44271) reduces to **just two
equations**:
```
    x_9770(A) = x_18274(B)      and      x_3183(A) = x_17728(B)
```
(atom 40782 is *implied* by these — proven in `test_40782.py`).

- `A` = the 22 control bits `BITS22` (drive x_9770, x_3183 only). **Fully enumerable: 2^22.**
- `B` = the other 233 control bits (drive x_18274, x_17728 only). 2^233 — the hard side.
- `x_18274 = x_6773/x_8821`, `x_17728 = x_17233/x_8821` (SHARED denominator x_8821).
- `x_8821` is **exactly linear** in the 233 bits; numerators are high-degree.
- best_partial_39019 sets ALL 255 control bits = 0.
- twist eqs: 1817 = 6033033*(x_9770-x_18274)+x_26977; 44271 = x_3183-x_17728;
  30378 = x_3183-x_9982-x_17728. (x_26977, x_9982 identically 0.)

## How to evaluate (the correct model)
`confluent_eval5.build5()` -> (A_atoms, kind, info, seq, bestval, ncyc). Build `seq`:
```python
order = json.load(open('eval_order.json'))['order']
defset = set(v for v in kind if kind[v] != 'const')
seq = [v for v in order if v in defset and v not in (9770,3183)]
seq += [v for v in (9770,3183) if v in defset]
seq += [v for v in defset if v not in set(order) and v not in (9770,3183)]
```
`make_forward(kind,info,seq,bestval)` -> Z solver `solve(list(bestval), setbits)`;
`make_forward(...,mod=P)` -> mod-P solver. forward_Z([]) violates exactly {1817,30378,40782,44271}.
The forward-eval satisfies every ORIENTED gate/load/div atom by construction for ANY bit set;
only the twist "check" atoms float — so it is a valid oracle for x_9770/x_3183/x_18274/x_17728.
NOTE: integer forward-eval is *lossy* (leaves a stale value when a division isn't exact) — use
the mod-P solver for any linearity/degree probing.

## Highest-EV next experiments
1. `runs/tab22_full.log` — full 2^22 (x_9770,x_3183) mod two 31-bit primes; saves
   tab22_9770_{p}.npy / tab22_3183_{p}.npy. When done: confirm B=0 fails; hash S and inspect
   structure (common factors, moduli). S then lets you INVERT the 22-side in O(1) (lookup).
2. Residue-pool identity: `extract_huge.py` -> huge_network.json (865 huge atoms; 512 simple
   loads bit*(x_B-HUGE)=s*x_C). Check whether x_9770(A) and x_18274(B) are combinations of the
   SAME HUGE residues => matching becomes combinatorial, not brute 2^233.
3. MITM/lattice via x_8821 (the linear coordinate on the 233 side) — see NOTEBOOK Session 6.

## Do NOT redo
- SAT/SMT (user directive: custom heuristics only; z3/cvc5 return unknown anyway).
- v4 evaluator / anything freezing x_18274 (fixed in v5).
- The lossy-eval "232-part slaved / rank 233" reduction (Session 5) — shown UNRELIABLE this session.
- Local search / greedy / pairs / triples from all-0 — all plateau (all-0 is the local min = 4 atoms).

## Git
Branch `claude/read-prompt-5t2raw`. Commit+push after meaningful experiments.
