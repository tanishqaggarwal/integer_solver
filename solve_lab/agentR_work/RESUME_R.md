# RESUME_R — agent R (automated reasoning against the REDUCED problem)

## 0. Score
- Baseline re-verified by me: `python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
  -> `satisfied 39026/39033 (7 failing)`, failing `[12231,12270,12350,14584,18673,22044,29125]`.
- **I have not beaten it.**

## 1. THE REDUCED PROBLEM IS NOW EXPLICIT (this is my main result)
Built from agent F's decode artifacts (read-only), all steps re-derived and checked here:

- `model.py` — p, offset K, law, target. Shift `X = x + K/3` removes the offset: the law becomes
  the plain chord law. Measured: **all 248 forced pin points and the target satisfy one and the
  same relation `Y^2 = X^3 + B`** with
  `B = 64019533680030876408443198762210829058751700634554282185987325820393598524794`
  (fitted from 2 points, then checked on all 248 + target: 0 exceptions).
- `group.py` — the fold is therefore an **abelian group law**: verified commutative and
  associative on 200 random triples, and `chordK` in the original coordinates agrees with the
  shifted chord law on 200 random pairs. **Consequence: the fold of a leaf subset does not depend
  on the tree shape at all.** F's remaining decode (56 slot pairs, 24 leaf-adjacent stages) is
  *not needed* to state or attack the problem.
- `ladder.py` -> `ladder.json` — the 256 conditional-pin leaves form **exactly one doubling chain
  of length 256**: 9 doubling-closed pieces (111,61,28,12,11,9,8,5,3 = 248 named) spliced by 8
  gaps of exactly one unknown point each. Checked: `L_i == 2^i * L_0` for sampled i.
- `order.py` — group order by Cornacchia (`4p = L^2 + 27M^2`), the unique candidate that kills
  the base point: **N = 115792089237316195423570985008687907852837564279074904382605163141518161494337**,
  a **256-bit prime** (`sympy.isprime`). Cyclic of prime order, so no Pohlig-Hellman.

**Therefore the reduced problem is exactly:**
> find the 256-bit integer `k` with `k * L_0 = T`, where bit *i* of `k` is the selector of
> ladder leaf `L_i` and `T = (30121525689829097248416773597728729849687459852468451992398421980273013515302,
> 44859544763832475231923253825569092119321525945631045653619508440821028887)` in shifted coords.

### F's requested validation — PASSED
Deliverable ON-set `{24601, 2081}` -> ladder indices **{72, 235}** -> `k = 2^72 + 2^235`.
`fold(k) != T`, exactly as F predicted the correct evaluator must report. The deliverable
*pins* the root wires to `T` while its own leaves fold to `2^72+2^235`; that gap is the 7 failing
equations.

## 2. Files
`model.py group.py ladder.py order.py ladder.json points_short.json` + `runs/` + `encodings/`.

## 3. Status / next
See LOG.md.
