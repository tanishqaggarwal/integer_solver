# Agent B — RESUME (independent re-parse track)

## VERIFIED
- `best/new_instance_partial_39026.json` -> 39026/39033, fails [12231,12270,12350,14584,18673,22044,29125]. CONFIRMED.
- Built my own parser+model, **verified exact** (`bverify_model.py`: 0/39033 mismatches on random ints):
    bparse3.py  -> parser keeping paren groups
    bmodel5.py  -> model5.pkl  (THE model; ~50 s to rebuild from source)
    bverify_model.py <seed>    -> proves the model reproduces every raw equation
    bscore.py <assign.json>    -> fast model score + residual decomposition
- **Every equation is `scalar * L^k = 0`, i.e. exactly `L = 0`**, L a linear form (|coef|<=80)
  over 38,133 gate atoms.  Kinds: sL 18478, pow2 8927, same 5821, plain 4256, pow4 783, same_pow2 768.
- **Only 5 factor shapes**: 3160 `x`, 9782 2-term-linear, 9067 3-term-linear, 16720 deg-2 2-term,
  512 deg-2 3-term.  39,241 distinct gate factors over 38,748 variables.
- p = 2^256-2^32-977 (secp256k1) is pinned by exactly ONE gate `x26064 - p` (f37650) and copied to a
  220-variable class.  3,707 gates have shape `y - q*p` with q occurring NOWHERE ELSE (free quotient).
- **ONLY FOUR big-constant pins exist in the whole instance**: f37650 (x26064=p),
  f39240 (x24453=97553848499418123410591666447050222001188385549510401465815187079080512838891),
  f1798 (x24468-x32989=91416...002), f729 (8863713*x18956 - x14257 = 11149...338).
- Max bipartite matching gate->output: 31,337 => 7,904 assertion gates, 7,411 free inputs.
- The 39,026 witness is DEGENERATE: 35,208 of 38,748 vars are 0, 2,555 are 1.
- **THE ENTIRE RESIDUAL = 4 divisibility conditions** (my own derivation, matches nothing inherited):
    (A) p | x28730   where x28730 = x4432 - x19964     [handle x9413, occurs in 1 gate]
    (B) p | x8731    (given x7075=1; x8731 = x3349 - x14865)   [handle x10903]
    (E) p | x9118    (x9118 = x32010 - x31861 = x31861 - x34310) [handle x1329]
    (D) 7376877*p | (x7068 - x2099)                    [handle x17325 via x642]
  x7075 = 1 - x21279.  Setting x21279=1 kills (B),(E) but perturbs x2099 via x25297=x9118*x21279.

## Best artifact
None above 39,026 yet.  Reference: /home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json

## Re-enter
    cd /home/user/integer_solver/solve_lab/agentB_work && python3 bmodel5.py && python3 bscore.py <json>

## Next experiment
Backward cone trace of x4432,x19964,x3349,x14865,x32010,x31861,x7068,x2099 to find free inputs that
move those residues mod p; then solve the 4 conditions jointly (CRT over the free handles).

## BREAKTHROUGH (agent B): exact circuit generator
`borient2.py excl_base.pkl orient3.pkl` + `ORIENT=orient3.pkl python3 beval.py <json> <out>`
recovers an ACYCLIC circuit orientation and forward-evaluates it:
  DEFINED = 34,184 gates,  FREE INPUTS = 4,557 vars,  ASSERTION gates = 2,767,
  excluded (allowed-nonzero) gates = 2,290 (the 2,283 dead branches of 2-factor atoms + the 7 residual).
Forward evaluation from the witness's OWN free-input values reproduces the witness EXACTLY:
  0 non-integral divisions, 0 differing vars, score 39,026, 0 violated assertions.
=> the whole instance is now a function of 4,557 integer knobs; all constraints are
   2,767 assertion gates + 7 residual divisibility conditions.

## Orientation refinement (agent B, later)
`borient7.py` + `bloop.py` iteratively forbid bad free-input promotions:
  DEFINED=34,253  FREE=4,488  ASSERTIONS=2,698  -> orient7.pkl
  ORIENT=orient7.pkl python3 beval.py <json> <out>  reproduces the witness exactly.
Cost of zeroing each residual congruence from the witness (ORIENT=orient7.pkl, bfix4.py):
  A p|x28730 via free x4432 : 39033-38995 = 38 failing
  B p|x8731  via free x3349 : 23 failing
  E p|x9118  via free x9118 : 24 failing
  D 7376877p|(x7068-x2099) via free x7068 : 27 failing
  all four   : 46 failing.   Baseline 7.
So each congruence is individually satisfiable; the cost is COLLATERAL in other gates.
Diagnostic tool: `ORIENT=orient7.pkl python3 bdiag.py` prints exactly which atoms break.
