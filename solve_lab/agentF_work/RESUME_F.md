# RESUME_F — agent F.  Self-contained handoff: rebuild everything from here.
No curve/group framing anywhere; integers, congruences and polynomials only.

--------------------------------------------------------------------------------------------------
## 0. SCORES (nothing here beats the shared baseline)
- Shared baseline **39,026 / 39,033**, re-verified by me:
  `python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
  -> `satisfied 39026/39033 (7 failing)`, failing `[12231,12270,12350,14584,18673,22044,29125]`.
- My own pipeline's best: **39,024** = `agentF_work/best_F_39024.json` (checker-verified, 9 failing).
  Also `best_F_39022.json` (39,022, only 2 nonzero atoms in the whole instance), `F_frame.json` (39,023).
- **No infeasibility is claimed.**  I held one two rounds ago and withdrew it; section 7 of LOG.md is
  marked INVALID.  Do not resurrect it.

--------------------------------------------------------------------------------------------------
## 1. WHAT THE INSTANCE IS  (all measured, all reproducible)
Every equation LHS is `(scalar)*S^k` or `(c1+c2)*S`, so each equation holds iff its **core S** is 0.
S is a left-nested spine of **atoms**; there are exactly **39,033 distinct atoms**.
30,001 atoms are *definitions* `x_out - f(inputs)` forming a **DAG**; forward evaluation from the
**8,747 free inputs** zeroes them all, leaving **9,032 residual atoms**.

**rank(M) = 39,033, dim ker(M) = 0** for the 39,033 x 39,033 equation-atom incidence matrix (525,982
nonzeros).  Three independent confirmations: a characteristic-free peeling certificate over Z
(`peel_order.npy`, re-verified by `peel_cert.py`; all pivots are +-1 or +-2, none divisible by any odd
prime, so it holds over Z and every field of char != 2), and Wiedemann over q=2^31-1 (minpoly degree
39,033, trailing coeff 268435456) and q=2147483629 (trailing coeff 11716781).
=> **any assignment satisfying all 39,033 equations must make all 39,033 atoms exactly zero.**
All-atoms-zero is EQUIVALENT to a full solve; a cancelling nonzero residual cannot exist.

Every check has the exact shape `A - B = M*p*h` with h a free integer, M a ~24-bit literal, and
**p = 115792089237316195423570985008687907853269984665640564039457584007908834671663** a 256-bit prime
literal of the file (it defines x26064).  So every check is a congruence mod M*p.
**The lift obstruction is at exactly ONE modulus, p, at level p^1 only:**
 - `modq.py`/`modm.py` fully solve the entire 39,033-equation system (0 nonzero atoms, 0 failing
   equations, ~3.5 s) modulo every prime tried from 3 to 2^255-19, every prime power tried
   (2^100, 3^80, 5^40, 7^25, 11^20, 13^20, 1009^8, 65537^4, 1000003^3, (2^31-1)^2), and the handle
   multipliers M themselves.  58 checkpoints in `modm_results/`.
 - mod p: 21 nonzero atoms / 124 failing.  mod p^2: 15 / 102.  Same code, same booleans.
 - `relaxed_pin.json`: relax the two unconditional constant pins and everything else closes exactly over
   Z -- the handles absorb every quotient.  So a mod-p solution lifts to Z immediately.
 - CRT of 20 independent 7-digit-prime solutions gives a solution modulo a **399-bit** composite Q with
   0 failing equations mod Q (`crt.py`, `crt_sols/`); its balanced representative scores 38,991.

--------------------------------------------------------------------------------------------------
## 2. THE DECODED STRUCTURE  — a 96-stage combination tree
**96 stage gates**, each carrying exactly 3 checks, each with its OWN distinct six-tuple of free inputs.
Root = gate **x15298** (leaf support 256).  Containment depth 6.  `tree96.json`.

**One uniform law, one universal constant.**  For every stage the three checks are satisfied exactly when
the output pair equals the chord-with-offset of the two input pairs:

    l = (b_y - a_y)/(b_x - a_x)
    out_x = l^2 - a_x - b_x - K
    out_y = l*(a_x - out_x) - a_y                            (all mod p)
    K = 97553848499418123410591666447050222001188385549510401465815187079080512838891

`stage_law2.py` searched all role partitions AND all coordinate orderings per stage, solving the 3x2
system once per partition and requiring two independent random draws to agree:
**72 of 72 stages with a full six-tuple obey it, and the offset is the SAME K in all 72 — zero
exceptions.**  Per-stage roles (which wires are inputs, which output, and the coordinate ordering) are in
`stage_roles.json`.  The remaining 24 stages have six-tuples with only 4 or 2 free inputs -- leaf-adjacent
stages whose missing inputs are literals hard-wired in the circuit; NOT yet resolved.
Stage x24533 was decoded by hand and its demanded output matched chordK(A,B) **digit for digit**.

**THE LAW IS INVERTIBLE — YES** (200/200 random triples, both directions, exact, O(1) each):

    B from (A,O):  l = (o_y + a_y)/(a_x - o_x);  b_x = l^2 - a_x - o_x - K;  b_y = a_y + l*(b_x - a_x)
    A from (B,O):  l = (o_y + b_y)/(b_x - o_x);  a_x = l^2 - b_x - o_x - K;  a_y = b_y + l*(a_x - b_x)

So the target inverts DOWN the tree as cheaply as leaves fold up, at every stage, not just the root.

**Wiring.**  Each stage input slot is a free variable w with a residual atom `((w - z) - handle)`; z
unfolds to a sum of gated `selector*value` terms.  `mux.py` (uses |Z| = 7,202 wires that are == 0 mod p
for EVERY assignment, to strip handles) decodes the 288 slot wires of the 72 full stages:
    116 -> 3-term gated mux ;  64 -> 1-term source ;
    108 -> not mux-fed: their atom shapes identify them as conditional LEAF pins `((X*(X-C))-(C*X))`,
           `((X*(X-C))-X)` or as stage CHECK atoms (output wires).
`mux_wiring.json` covers **47 of 72** stages.  **56 stages still have an undecoded slot pair and NOTHING
was filled by guess.**

**The 3-way mux is the three non-zero quadrants of two OR-groups** -- verified:
    x34606 = x7715*x23597 = a*(1-b),  x5647 = x34554*x19271 = b*(1-a),  x15298 = x7715*x34554 = a*b
so **exactly one branch is live** whenever the subtree has any live leaf.

--------------------------------------------------------------------------------------------------
## 3. PRIORITY 1 ANSWERED — the reachable space is 2^256-1, and I was WRONG to suggest otherwise
Because exactly one mux branch is live per slot, EVERY boolean assignment yields a well-defined fold.
The count then obeys, for a stage with child slots S1,S2 (fold when both live, pass-through when one is):
    N(stage) = N(S1)*N(S2) + N(S1) + N(S2),   N(leaf) = 1
which solves to **N(n) = 2^n - 1** exactly for a subtree with n leaves (N(2)=3, N(3)=7, ...).
So the reachable configuration count at the root is **2^256 - 1 non-empty leaf subsets**.
**This CORRECTS my previous round.**  I warned that the phrase "choose a subset of 256 leaves" might be
wrong because the gating could restrict reachable patterns.  It does not restrict them -- the quadrant
structure makes every non-empty subset well-defined.  The original phrasing was right; my caution was not.
(Caveat, stated: this is a model from 47/72 wired stages plus the verified quadrant structure, not an
exhaustive check.  What remains genuinely open is what happens when two leaves in the SAME OR-group are
both ON -- both pins fire on different wires and the slot may then see a sum rather than a single value.
That case is not covered by the count and must be settled by the completed decode.)

**Consequence: enumeration is NOT the attack.  Inversion is.**

--------------------------------------------------------------------------------------------------
## 4. PRIORITY 4 DATA — leaf-support profile of all 96 stages (`stage_profile.json`)
Root split: inA = (x12186,x16742) <- sources (x23927,x19083), **leaf support 178**;
            inB = (x14853,x24908) <- sources (x1308,x17601),  **leaf support 78**; out = (x22162,x30213).
Only **three** stages exceed 24 leaves: **256 (x15298), 88 (x30973), 50 (x24533)**.
**93 of 96 stages have leaf support <= 24**, and 66 lie in the window 2..24.
Memory ceiling on this box: 2^27 entries x 64 B = 8.6 GB, so realistically **2^24-2^25 entries (1-2 GB)**.
A meet-in-the-middle at the root needs 2^78 on the smaller side -- hopeless.  But the window is populated:
**invert the target down through the 78-side chain to a node of support <= 24 and enumerate there.**
The chain from the root to those nodes runs through x30973 (88) and x24533 (50), so the inversion has to
pass two large nodes; each pass is O(1) per candidate by the formulas in section 2.

--------------------------------------------------------------------------------------------------
## 5. WHAT IS NOT DONE (believe nothing from these until they exist)
- The **fold evaluator is not built and not validated**, and **no subset search has been run**.
- 56 of 72 stages have an undecoded slot pair; the 24 leaf-adjacent stages' literal inputs are unresolved.
- The same-OR-group double-leaf case (section 3 caveat) is unsettled.
- Agent E's 7 residue classes were NOT compared with my channels: with 56 stages incomplete a partial
  match would be spurious.  Do it after the decode completes.

## 6. NEXT, IN ORDER
1. Finish the decode: the 56 stages and the 24 leaf-adjacent ones.  Settle the same-group double-leaf case.
2. Build the evaluator over `tree96.json` + `stage_roles.json` + `mux_wiring.json`.
   **VALIDATE**: reproduce the deliverable's ON-set {24601, 2081}; the prediction that must hold is
   root = fold of those two leaves and **NOT** the target.  If it equals the target, the evaluator is wrong.
3. Deep meet-in-the-middle: invert the target down the 78-side to a node of support <= 24, enumerate that
   side forward, match in a hash table of <= 2^24 entries.
4. Then compare channels with agent E's residue classes.

## 6b. HANDOFF VERIFIED (run on the last session, not asserted)
All 42 files and both directories named below exist.  `python3 parse3.py; python3 circ4.py;
python3 sched.py; python3 supp.py` runs clean from cold.  `python3 fwd.py` reproduces
"3 nonzero residual atoms, 28 failing equations => 39005"; `python3 peel_cert.py` reproduces
"certificate verified: True, 39033 of 39033, rank(M)=39033, dim ker(M)=0"; `checker.py
best_F_39024.json` reproduces 39024/39033.  The ONE thing this document describes that does not exist is
the fold evaluator, and it is documented as not existing (section 5).
Note on a withdrawn fleet-wide criterion: agent H's "rank > deficit" test is retracted.  I never used it --
my frame pricing is integer reachability via `intsolve.solve_int` (column HNF), and rational rank appears
in my log only as the quantity that fails to predict.  Nothing of mine needs re-checking on that account.

## 7. FILES
Code: `parse*.py circ*.py sched.py supp.py fwd.py full.py jac.py intsolve.py lin.py frame.py gs2.py
modq.py modm.py crt.py buildM.py peel.py peel_cert.py wiedemann.py cfg_rigid2.py stage_law2.py mux.py`
Data: `tree96.json stage_roles.json mux_wiring.json stage_profile.json sweep_ii.json sweep_i.json
M.npz peel_order.npy modm_results/ crt_sols/ supp.pkl circ4.pkl sched.pkl`
Rebuild pickles (~1 min): `python3 parse3.py; python3 circ4.py; python3 sched.py; python3 supp.py`
Full narrative with every measurement and every retraction: `LOG.md` (sections 1-32).
