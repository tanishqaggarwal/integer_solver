# The obstruction as three congruences (CRT reduction)

Best verified: `best_agentA_39022.json` — 39022/39033, 11 fails
`[2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125]`.

## All 11 fails = {G1=0, G2=0}
Every failing equation is a fixed linear combination of two gap conditions (the
gate-definition atoms hold automatically under forward evaluation). E.g. eq 22044
(atom 42669) = `-24·G1 - 29·G2` once the gate atoms A,B,C vanish.

```
G1 = x_7068 - x_2099 - 7376877 * x_642 ,   x_642   = x_28599 * x_17325
G2 = x_4432 - x_19964 - x_28730       ,   x_28730 = x_17499 * x_9413
```
- `x_28599`, `x_17499` are gate outputs **rigidly = p** (copy chains from the
  hardwired constant gate `x_26064 = p`; `x_28599 = x_9325 = … = p`). Not movable.
  ⇒ `x_642 = p·x_17325`, `x_28730 = p·x_9413` (p-granular slacks).
- `x_2099 = x_6418` pinned to constant `M1` by atom 3277 = `x_2081·(x_6418-M1)`
  (active since selector a=x_2081=1). `x_19964 = x_12553` pinned to `M2` by atom 3279.
- `x_7068, x_4432, x_17325, x_9413, x_6418, x_12553` are all **free inputs**.

## The absorbers are side-effect-free (verified)
Setting `x_17325`, `x_9413` to ANY value (tested up to 7·10⁴⁰) keeps exactly the
same 11 fails — they appear only in G1/G2's atoms + the master combiner. So once
the congruences below hold, the absorbers close G1,G2 with no collateral.

## The reduction
Solving G1=G2=0 over ℤ is EXACTLY:
```
(A)  x_7068 ≡ x_2099  (mod p)            [= M1 mod p]
(B)  x_7068 ≡ x_2099  (mod 7376877)      [= M1 mod 7376877;  7376877 = 3²·819653]
(C)  x_4432 ≡ x_19964 (mod p)            [= M2 mod p]
```
then set `x_17325 = (x_7068-x_2099)/(7376877·p)`, `x_9413 = (x_4432-x_19964)/p`
(both exact integers once A,B,C hold).

G2 needs only a mod-p alignment (its slack step is p, absorber x_9413 takes any
carry). G1 additionally needs the small-modulus condition (B) because its slack
step is 7376877·p. M1, M2 are hardcoded constants, so x_7068 must be *moved* to
≡ M1 (mod 7376877·p) and x_4432 to ≡ M2 (mod p) with their gate-subsystems kept
consistent. Current gaps: (x_7068-x_2099) mod 7376877 = 3228258;
mod p ≠ 0; (x_4432-x_19964) mod p ≠ 0.

## REFINED (the pins have p-absorbers → obstruction collapses to mod p only)
The message pins are NOT rigid. Their atoms carry absorber terms, and the absorber
multiplier wires equal p:
- atom 3277 (satisfied eq): `x_2081·(x_6418-M1) - 15804267·x_26777 = 0`,
  `x_26777 = x_38744·x_3387` with `x_38744 = p`. ⇒ `x_2099 = M1 + 15804267·p·x_3387`.
- atom 3279: `x_2081·(x_12553-M2) - x_13458 = 0`, `x_13458 = x_22972·x_5081`,
  `x_22972 = p`. ⇒ `x_19964 = M2 + p·x_5081`.
VERIFIED side-effect-free: moving x_3387 (with x_6418 following the pin) or x_5081
(with x_12553 following) keeps exactly the same 11 fails — zero collateral. So there
are FOUR clean absorbers: x_3387, x_17325 (G1 side), x_5081, x_9413 (G2 side).

All four absorber shifts are multiples of p. Hence they can hit every residue class
EXCEPT mod p. Therefore:
- For any modulus m coprime to p (3, 9, 819653, …): the gap condition is absorbable
  (e.g. x_2099 mod 819653 is fully tunable via x_3387, p a unit). NOT binding.
- **The only binding obstruction is mod p:**
  ```
  (A)  x_7068 ≡ M1 (mod p)          (C)  x_4432 ≡ M2 (mod p)
  ```
  (M1 = x_2099 mod p, M2 = x_19964 mod p — both fixed mod p.) Once A,C hold, set the
  four absorbers (with mild 3 | (x_7068-M1)/p via Bezout, gcd(15804267,7376877)=3) to
  zero G1,G2 with no side effects.

Current binding residuals in best_agentA:
  (x_7068-M1) mod p = 61705020361863629770768910187978745858728889529652486596432934143473517757811
  (x_4432-M2) mod p = 33310166114805471624282140578459083391052142224394967852279417483154815501175

So the whole 39033 problem = "can channel (1,0) be solved mod p with the response wires
x_7068, x_4432 matching the message residues M1, M2 mod p, keeping the other 39022 eqs
satisfied?" A global GF(p) linear-feasibility (Newton from agentA). Feasible ⇒ lift +
absorbers ⇒ solved. Infeasible ⇒ channel (1,0) is ℤ-infeasible ⇒ correct MUX msg differs.

## Consequence for method
An ℤ solution ⇒ a solution mod every m. So channel (1,0) feasibility can be tested
cheaply mod small moduli:
- mod 3:  7376877≡0 ⇒ G1≡x_7068-x_2099; already ≡ (both 2 mod 3).
- mod 9:  G1≡x_7068-x_2099, currently off by 3; G2 mod 9 absorbs via x_9413.
- mod 819653: G1≡x_7068-x_2099 (7376877≡0), need x_7068≡M1; G2 absorbs.
- mod p:  slacks vanish, need x_7068≡M1 and x_4432≡M2.
If channel (1,0) is infeasible mod any small m → it is ℤ-infeasible and the correct
MUX channel differs. Channels (1,1)/(0,1)/(0,0) give 30/25/129 naive fails; (1,0)
is the unique best basepoint.
