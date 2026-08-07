# Agent B — RESUME (independent re-parse track).  FINAL.

## Best verified score: 39,026 / 39,033  (no improvement over the lab's deliverable)
- Reference: `/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json`
  re-verified by me with `solve_lab/checker.py` -> `satisfied 39026/39033`,
  failing `[12231,12270,12350,14584,18673,22044,29125]`.
- My own equivalent state, produced entirely by MY generator (not copied):
  `solve_lab/agentB_work/out/fwd7b.json` — checker-verified 39,026/39,033, same 7 failing lines.

## Pipeline (all derivable from source; every .pkl is a regenerable cache)
    cd /home/user/integer_solver/solve_lab/agentB_work
    python3 bmodel5.py                 # ~50 s -> model5.pkl        (parser: bparse3.py)
    python3 bverify_model.py 4242      # ~35 s -> "TOTAL mismatches: 0 / 39033"  <-- PROOF the model is exact
    python3 bscore.py <assign.json>    # fast model score + residual decomposition
    # rebuild the circuit orientation:
    python3 -c "import sys,pickle;sys.path.insert(0,'.');import os;os.environ['ORIENT']='orient7.pkl';\
import beval as E;v0=E.load('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json');\
pickle.dump([f for f in range(len(E.facs)) if E.fval(f,v0)],open('excl_base.pkl','wb'))"
    rm -f forbid.pkl && python3 bloop.py       # ~2 min -> orient7.pkl
    ORIENT=orient7.pkl python3 beval.py <assign.json> <out.json>   # forward-evaluate the circuit
    ORIENT=orient7.pkl python3 bfix5.py        # drive the four residual congruences
    ORIENT=orient7.pkl python3 bdiag.py        # print exactly which atoms break
    ORIENT=orient7.pkl python3 bscan2.py       # (written, not finished) derivative scan of all free inputs

## ESTABLISHED — my parser shares NO code with the prior lab
1. `bverify_model.py` proves my decomposition reproduces every one of the 39,033 raw equations
   EXACTLY on random integer assignments: 0 mismatches, 2 independent seeds.
2. **Every equation is `scalar * L^k = 0`, i.e. exactly `L = 0`**, with L a linear form
   (|coef| <= 80) over 38,133 gate atoms.
   Kinds: sL 18478, pow2 8927, same 5821, plain 4256, pow4 783, same_pow2 768.  No exceptions.
3. 39,241 distinct gate factors over 38,748 variables.  Only FIVE factor shapes:
   3160 `x`; 9782 2-term linear; 9067 3-term linear; 16720 deg-2 2-term; 512 deg-2 3-term.
   48 anonymised gate shapes: mul, square, add/sub, copy, const-mul, boolean, const-pin,
   plus conditional pins `b*(x-K) = c*handle`.
4. p = 2^256-2^32-977 (secp256k1) is pinned by exactly ONE gate `x26064 - p` (f37650) and copied
   to a 220-variable class.  3,707 gates have shape `y - q*p` where q occurs in NO other gate.
5. Exactly FOUR big-constant pins in the whole instance:
   f37650 `x26064 = p`;
   f39240 `x24453 = 97553848499418123410591666447050222001188385549510401465815187079080512838891`;
   f1798  `x24468 - x32989 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002`;
   f729   `8863713*x18956 - x14257 = 1114942656963403660822546820446916783439088877768247923308647546252105232931473698035897478439338`.
6. PARSE KEY: the generator wraps the HEAD gate of every packing chain in an EXTRA paren group.
   That is the reliable gate boundary (`bparse3.flat_pack`).  Without it the decomposition
   over- or under-splits (94k / 58k spurious atoms instead of 38k real ones).
7. **CIRCUIT GENERATOR** (`borient7.py`/`bloop.py` -> orient7.pkl): acyclic orientation with
   34,254 defining gates, 4,487 free inputs, 2,697 assertion gates.  Forward-evaluating from the
   witness's OWN free-input values reproduces the witness EXACTLY: 0 non-integral divisions,
   0 differing variables, score 39,026, 0 violated assertions.
   => the entire instance is a function of 4,487 integer knobs.
8. The 39,026 witness is DEGENERATE: 35,208 of 38,748 variables are 0, 2,555 are 1.
9. **THE WHOLE RESIDUAL = 4 divisibility conditions** (my own derivation):
     A  p | x28730,  x28730 = x4432 - x19964,  x19964 = x20492 = x36065 = x12553  [handle x9413]
     B  p | x8731,   x8731  = x3349 - x14865                                       [handle x10903 via x31864]
     E  p | x9118                                                                  [handle x1329  via x29854]
     D  7376877*p | (x7068 - x2099),  x2099 = x6418                                [handle x17325 via x642]
   Free knobs: A <- x4432 or x12553; B <- x3349; E <- x9118; D <- x7068 or x6418.
10. MEASURED COSTS from the witness (`ORIENT=orient7.pkl python3 bfix5.py`):
      B+E together: BOTH congruences zeroed with **ZERO collateral** -> only 3 nonzero atoms
                    remain (A and D), score 39,023 (10 failing).
      A via x12553 : 39,009.   D via x6418 : 39,008.
      A+B+D+E via x12553/x6418 : 39,019 (14 failing) — all four congruences zeroed simultaneously.
      A+B+D+E via x4432/x7068  : 39,001, and the ONLY nonzero atoms are the three collateral ones.
    The three collateral congruences that block a full solve:
      G1  p | (x24548 - x25442)        [f8715: x7927 = 9367949*(x24548-x25442); f8714: x7927 = x11052*p]
      G2  8481759 | (x15324 - x37254)  [f22518: x15324 = x4432 + x10250; f11925 handle x37413]
      G3  p | (x2964 - x26756)         [f8713: x579 = x2964 - x26756; f20043: x579 = x19569*p]
    Exact derivatives measured (`bmeas.py`):
      d(x7927)/d(x4432) = -9367949 (mod p);  d(G2)/d(x4432) = +1;  d(x579)/d(x7068) = -1 (mod p).
      Shifting x4432 by p moves ONLY G2.  Shifting x7068 by 7376877p, x3349 by p, or x9118 by p
      moves NOTHING.  So the congruence lattice on the knobs is fully explicit.

## SINGLE HIGHEST-VALUE NEXT EXPERIMENT
Solve {A,B,E,D,G1,G2,G3} as ONE CRT / linear-Diophantine system in the free knobs
{x4432, x12553, x3349, x9118, x7068, x6418, x2964, x10250, x24548, x25442, ...}.
`bscan2.py` is already written and fills in the full 7 x 4,487 derivative matrix
(one forward evaluation per free input, ~0.2 s each, ~20 min total).
B and E already cost NOTHING, so only A and D need a collateral-free driver; the search
space is small and fully characterised above.
