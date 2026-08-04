# RESUME — read me first

## ⚠️ CURRENT (re-randomized) INSTANCE — NOT solved; identified as a secp256k1 GF(p) trapdoor
The EQUATIONS.txt in the repo is a NEW re-randomized instance (39,033 eqs). Full analysis in
**`NEW_INSTANCE_STATUS.md`** (read it). Best verified partial: **39,007 / 39,033**
(`best/new_instance_partial_39013.json`).

Definitive findings (exhaustive):
- Gate DAG fully ACYCLIC; forward-eval from free inputs satisfies all wiring automatically.
- The ONE large identity wire (220 vars) is PINNED to **p = 2^256-2^32-977 = secp256k1 field
  prime** (x_26064=p, appears 13x); it is the twist multiplier. No free wire exists (unlike the
  solved instance, whose wire was FREE — that was the whole solve).
- The remaining obstruction is a **256-bit boolean codeword message** (disjoint 178+78-bit cones
  of control bits x_7715,x_34554). x_9274=OR(controls)=1 is FORCED, so activation is mandatory.
- Each set bit triggers a huge additive load; the GF(p) load matrix over the 256 bits is FULL
  RANK (no bit self-cancellation). Data can absorb mod p but the ℤ-lift imposes p-divisibility
  carries -> every wrong-message data solve is ℤ-inconsistent (SNF pivot=p); iterative repair
  diverges (27->300+). Inhomogeneous GF(p) message solve inconsistent (rank 3). No small vinegar
  linearization (21,922 vars bilinear). => structured GF(p) codeword/MQ; needs the setter secret.
- Key tools: build_twist.py (activate+route MUX), newton2.py (simultaneous absorber solve),
  p_message.py (GF(p) message solve), scc.py, localize.py, scan_bits.py, p in huge_consts.json.
- Untried heavy attack: lattice/LLL on the mod-p codeword (standard, likely resisted by design).

---
# RESUME — read me first

## ✅ SOLVED — all 39,031 / 39,031 equations satisfied exactly in ℤ
Solution: `best/SOLUTION.json` (also repo root `SOLUTION.json`), 38,748 vars, 2,954 nonzero.
Verify: `python3 checker.py best/SOLUTION.json` → `RESULT: OK`.
Reproduce: `python3 build_solution.py`. Full writeup: **`SOLVED.md`**.

The twist reduced to activating two product-slacks (H: x_9982=x_12518·x_9897=−G;
F: x_26977=x_20510·x_31302=F0). The hubs x_12518/x_20510 sit in a 220-var identity wire
that is *exactly quiet*, so it is a free parameter; set wire=sign·1 and the two rare
partners x_9897=−G, x_31302=F0. Direct algebraic construction, no search. The earlier
forward-eval framing (control bits, div-wire, x_12779≥2, dirty bits, high-dim co-activation)
was an artifact of the orientation — none of it was needed.

---
## (historical) prior best partial
Best verified partial was **39,019 / 39,031** (`best/best_partial_39019.json`).

## SESSION 7 — SLACK-ACTIVE SOLVER BUILT; obstruction reduced to R=0 (read NOTEBOOK Session 7 tail)
The slack-active evaluator EXISTS now: `slack_active.py` (freeze x_24026:=x_18274-x_35186,
x_27116:=x_17728-x_1642 with x_12779=1 via a single 22-side bit e.g. 1858). It makes BOTH twist
halves hold by construction — the state plain forward-eval cannot represent. Activating the
slack ripples into ~18 verifier CHECK atoms; SA-with-square-roots (`slack_sa.py`, replaces the
deg-4 squares a40782/a39550 by their deg-2 roots Q=0 via `check_square.try_sqrt`) drives the
frustrated core 18 -> 6. Run the 4-way fleet: `python3 slack_sa.py <activator> <seed> <out.json>`
with activators in {1858,26947,27512,30104,5443,...}.
CRISP OBSTRUCTION: for verifier square a40782, satisfying it AND a1817 reduces to R=0 where
R = 28*x_10783 + (ripple terms), x_10783=x_16644*x_17301, all fixed by the RIGID 3183-slack
(a44271: x_3183=x_17728, so x_27116=x_17728-x_1642 is pinned). The continuous knobs x_24026 and
the FREE var x_31302 (df=None) CANNOT change R (Q40782 slope in x_24026 is 0 once a1817 held).
So the witness = a DISCRETE 233/22-bit choice whose rigid-slack ripple self-annihilates in every
verifier square (a knapsack). The div-wire escape (x_8821=x_17810*x_27292 in {-2,-1,0,1}) lets
x_18274/x_17728 leave their g2/h2 lattice but only onto (base/2)*Z, still coprime => degenerate.
NEXT: keep the slack-active SA fleet running; or attack R_i=0 across the 530 squares as a system
(linearize ripple monomials); or find the setter's 233-bit knapsack solution (LLL blocked by
numerator nonlinearity — 7/50 linear). NOTE: single 22-side bits give x_12779=1 (not 2).

## TRAPDOOR MECHANISM — fully reverse-engineered (Session 6, read NOTEBOOK Session 6 tail)
The obstruction (atoms 1817,30378,40782,44271) is the twist x_9770=x_18274 & x_3183=x_17728.
KEY: the confluent forward-eval QUANTIZES both sides to COPRIME units and ZEROS the slack
products, so it can NEVER represent the (feasible) witness — this is why every forward-eval
search (SA/mitm/greedy/pairs/enum/local) plateaued. Specifically:
- Under forward-eval: x_9770=m*g, x_3183=m'*h, x_18274=m2*g2, x_17728=m2'*h2 (g=119182..,
  g2=91416..; gcd(g,g2)=1, gcd(h,h2)=2). Rigid twist => degenerate 0 only. (codewords.py, quant_structure.py)
- BUT the wire DEFS carry product slacks: x_9770 = x_35186(=m*g) + x_3368, x_3368=x_12779*x_24026;
  x_3183 = x_1642(=m'*h) + x_10466, x_10466=x_12779*x_27116. Both gated by x_12779=x_23380*x_36336.
- forward-eval sets x_12779=0 (slacks off) -> quantization. The WITNESS activates x_12779 (22-side
  bit pairs give x_12779=2) AND x_24026/x_27116 (deeper, via x_38215) so
      x_9770 = m*g + x_12779*x_24026 = x_18274 = m2*g2   (bridges the coprime gap).
- So the TRUE solve = search WITH the slacks active. With slacks on, x_9770 is NOT limited to 27
  values and CAN equal x_18274; the decoupling (x_9770<-22 only) is a slacks-OFF artifact.

NEXT-STEP for a solver: build an evaluator/search that DRIVES x_12779, x_24026, x_27116 nonzero
(find their activating bit cascades: x_12779<-{1858,2795,5443,10652,19520,26947,27512,30104,...},
x_24026<-x_38215<-...), then solve m*g + x_12779*x_24026 = x_18274(B) (coupled product match).
Do NOT rely on the all-0 forward-eval regime — it structurally excludes the witness.

## Earlier (still true) reduction
- `A` = the 22 control bits `BITS22`; `B` = the other 233 bits.
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

## Exhausted this session (do NOT redo)
- SAT/SMT (user directive: custom heuristics only; z3/cvc5 return unknown anyway).
- v4 evaluator / anything freezing x_18274 (fixed in v5).
- Linear algebra / lattice: `linalg255.py` (CORRECT, over all 255 bits, mod-P) has RANK 255/255
  and forces ALL bits = 0. The witness (!= all-0) is OUTSIDE the linear neighborhood of all-0, so
  linear/lattice/subset-sum provably cannot reach it. Supersedes the Session-5 "slaved-B" claim.
- B=0: ruled out (full 2^22 scan, 0 matches).
- Modulus (gcd residues=1), residue-lattice relation (none), slack vars (x_26977/x_9982 rigid).
- Local search / greedy / SA / pairs / triples from all-0 — all plateau (all-0 is the local min).

## The ONE remaining avenue (unimplemented)
A custom NONLINEAR solver / backward circuit-inversion: pin x_18274=N1, x_17728=N2 for a chosen
(N1,N2) from the 2^22 table S, and propagate/search backward through the 233-side acyclic circuit
(residue-load selects + product/sum gates) to determine the bits. Big build, uncertain (z3 failed
the analogous forward CSP). This is the only path not proven dead — everything else is exhausted.
The instance is a genuine obfuscated-circuit trapdoor; a full witness likely needs the setter's
secret or a cryptanalytic break of the specific 233-side residue circuit.

## Git
Branch `claude/read-prompt-5t2raw`. Commit+push after meaningful experiments.
