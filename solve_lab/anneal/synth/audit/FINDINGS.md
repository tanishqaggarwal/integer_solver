# Adversarial audit -- FINDINGS

Role: try to BREAK faithfulness (a spurious ground state) and soundness (a false
accept). **Headline: no spurious ground state and no false accept found.** Five
surviving concerns, all LOW. Scripts in synth/audit/ (a2_modmul_multirhs.py,
a3_gadgets.py, a4b_ladder_feasible.py, a5_andweight.py, a6_reduction.py).

## Verified the verifier
Cross-checked verify.py's zero_states() against a true 2^n brute force over the final
Hamiltonian: p=3 n=21 bf_full=28 == zs_full=28. The enumerator returns exactly the
true zero-energy set; used as the instrument below.

## Attack 1 -- spurious ground states: NONE
- 1a. Modmul in the forms the ladder ACTUALLY uses (multi-word RHS + squaring
  lam^2==x3+x1+x2, lam*mm==y3+y1 -- untested by shipped verify.py): exhaustive over
  inputs AND full ground-state enumeration, every reduction/mode/const/square:
  GRAND TOTAL FAILURES 0, |E=0|==|truth|.
- 1b. Ladder gadgets (full enumeration): one-hot sequential counter zero-energy set =
  exactly the D one-hot vectors (no all-zero, no multi-hot); not_equal(W,c) =
  exactly {0..2^s-1}\{c}.
- 1c. DECISIVE ladder test (a4b): clamped digit inputs, asked whether ANY zero-energy
  completion exists (demo_win2 only replays the canonical witness). Feasible <=> true
  solution, exactly; the 'degen' rows (chain hits x1==x2) are INFEASIBLE -- the d!=0
  gadget blocks the teleport at the level of ALL completions. Strictly stronger than
  the shipped demos. TOTAL BREAKS: 0.
- 1e. L0X matched |E=0|==|truth| for p=3,5,7,13: quotient-word width neither too loose
  (spurious q) nor too tight (dropped states).

## Attack 2 -- AND-weight sufficiency
Faithful for any W>=1: every square >=0 and every Rosenberg AND penalty >=0 (=0 iff
z=ij), so any gate-violating state has E>=W>=1>0. W_and=1 gives the same E=0 set as
default. The local bound is COMPLETE because no AND output ever feeds an AND input
(degree <=2 monomials), so and_weight_ok() holds everywhere. No counterexample.

## Attack 3 -- reduction, independently recomputed from raw pinrec.json + huge_consts.json
Independent cubic fit: all 256 pairs on curve; matches core a2/a4/a6; target on curve;
A==0; B matches; doubling chain heads=1 len=256, P_i==2^i G for all i; n prime, n*G==O,
n*T==O. Everything the encoder consumes is reproducible bit-for-bit; chain exactly 2^i G,
no off-by-one.

## Attack 4 -- degenerate EC cases
Real instance: no P_i/T is O; none has y=0 (odd prime order => no 2-torsion); 256 chain
x-coords distinct; 2^255<n. Addend-level infinity/repeat/2-torsion cannot occur.
Accumulation-time x1==x2 possible but rendered INFEASIBLE by d!=0 (proven in 1c).

## Attack 5 -- end-to-end
No false accept: checker.py REJECTS all-zeros (11684 failing) and best_agentA_39022
(11 failing); demands all 39,033 equations exactly 0 in Z. Planted recovery works
(bits=16,24). Multiplicity concrete: bits=24 seed=3 (n<2^24) has TWO solutions k and
k+n, both <2^24, both mapping to T.

## Surviving concerns (all LOW)
1. (LOW-MOD) Solution multiplicity real and was under-stated: k and k+n both valid when
   dlog < 2^256-n ~ 2^128. No soundness impact (both true solutions); the "exactly ONE"
   wording in solve.py/demo.py/README was corrected to "1 or 2 (k and k+n)".
2. (LOW) Completeness gap from d!=0: a true k whose ladder hits a degenerate case has no
   E=0 state, so "no E=0 found" does not prove UNSAT. Bounded ~2^-247 (ENCODING.md).
3. (LOW) Shipped harnesses replay only the canonical witness; closed independently (1c).
4. (LOW) build_modmul multi-word-RHS/squaring were unverified by shipped verify.py
   though in fact faithful (1a); recommend adding those cases.
5. (INFO) Exhaustive proofs are toy-scale; the 256-bit claim rests on the compositional
   argument (each gadget uniquely pins its output mod p) being scale-invariant, which it
   structurally is (bit-widths/quotient ranges formula-derived).

## Bottom line
Could not produce a spurious ground state or a false accept. Modmul faithful in every
form the ladder uses; both ladder gadgets faithful; the full ladder admits an E=0
completion IFF the digit-vector is a true solution (degenerate collisions provably
blocked at the all-completions level); AND weight sufficient for any W>=1 with a complete
local bound; reduction independently reproducible; checker.py rejects non-solutions. The
certificate survives, subject to the honestly-scoped multiplicity (#1) and completeness
(#2) caveats.
