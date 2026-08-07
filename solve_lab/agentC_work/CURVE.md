# THE CURVE — read this before testing any point against `y^2 = x^3 + 7`

> **This lab has produced THREE false negatives from the same mistake: testing the instance's
> points against the SHORT form `y^2 = x^3 + 7` when the circuit's curve has a nonzero `a2`
> (the `x^2` coefficient). One of those false negatives stood for many sessions
> (`s10/curve.py`, "secp256k1 hypothesis refuted", Part VII). It is wrong. Do not repeat it.**

Regenerate every number below from `EQUATIONS.txt` alone:
`python3 agentC_work/CERT_second_door.py` (section 3 fits the curve; section 5 checks the order).

---

## 1. The form the circuit actually uses

The root point-addition residual measured off the instance is

```
A = (x2 - x1)^2 * (x3 + x1 + x2 + a2) - (y2 - y1)^2
B = (y3 + y1)*(x2 - x1) - (x1 - x3)*(y2 - y1)
```

with `x1 = x_12186, y1 = x_16742, x2 = x_14853, y2 = x_24908, x3 = x_22162, y3 = x_30213`.
These match `x_35389` and `x_6671` **digit for digit** over random probes.

The `a2` inside `A` is **not** an anomaly, an "extra constant K", or evidence against the
elliptic-curve reading. For the general Weierstrass form `y^2 = x^3 + a2 x^2 + a4 x + a6`
the chord-and-tangent law is exactly `x3 = lambda^2 - a2 - x1 - x2`. Setting `a2 = 0` and then
observing a leftover constant is simply reading the general law through a short-form lens.

## 2. The exact constants (p = 2^256 - 2^32 - 977)

```
a2 = 97553848499418123410591666447050222001188385549510401465815187079080512838891
a4 = 114170008767671698752186727197936107864370654164657728518655355473804451402762
a6 = 77755683306591771556999954628254672912734268662742093169295805431582354953490
```

Fitted from three pinned leaf constants (3 linear equations, 3 unknowns) and then verified
on the other 253: **all 256 free leaf points lie on this one curve.** Discriminant != 0.

Machine-readable copy: `agentC_work/curve.json` (keys `KA` = a2, `a4`, `a6`, `P1`, `P2`, `Q`).

## 3. Short form, and why it is secp256k1

Substitute `X = x + a2/3`:

```
A_short = a4 - a2^2/3 = 0        <-- EXACTLY zero, so j = 0
B_short = 64019533680030876408443198762210829058751700634554282185987325820393598524794
```

so the curve is `y^2 = X^3 + B_short`.

**Sixth-power witness.** `B_short / 7` is a sixth power mod p:

```
pow(B_short * pow(7, p-2, p) % p, (p-1)//6, p) == 1     -> True
```

Hence `B_short = 7 * u^6` for some `u` in F_p, and the map `x -> u^2 x, y -> u^3 y` is an
**F_p-isomorphism onto secp256k1**. It is the TRIVIAL sextic twist.

**Group order, pinned exactly.**

```
n = 115792089237316195423570985008687907852837564279074904382605163141518161494337
  = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
```
which is the secp256k1 group order. Verified: `n` is prime (Miller-Rabin, 40 rounds) and
`[n]G = O` for the chain generator `G`. Since `n` is prime and `G != O`, `ord(G) = n` exactly.

> A *nontrivial* sextic twist of a `j = 0` curve carries one of the other five CM orders
> `p + 1 - t` with `t` in `{+-L, +-(L+9M)/2, +-(L-9M)/2}`, `4p = L^2 + 27M^2`. Same-`b`-different
> is not the test; **same order is**. So "different `b`, therefore a twist" is the wrong
> inference — a different `b` with the *same* order is an isomorphic copy.

## 4. The test that produced the false negatives

`Q = (K2 mod p, K1 mod p)` where
`K1 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626`
`K2 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002`

| test | result |
|---|---|
| `Q` on `y^2 = x^3 + 7 (mod p)` | **False**  <- the false negative |
| `Q` on `y^2 = x^3 + a2 x^2 + a4 x + a6` | **True** |
| all 256 leaf points on `y^2 = x^3 + 7` | 0 / 256 |
| all 256 leaf points on the real form | **256 / 256** |

## 5. Consequences that depend on this being right

* The 256 free selector bits carry `P_i = 2^i G` (verified for **all** i = 0..255), so the
  selector tree accumulates `k*G` with `k` the 256-bit integer of the bits.
* A full solve needs either `P1 + P2 = Q` (the discrete log of `Q`, on a group F_p-isomorphic
  to secp256k1, order the prime `n`) or the degeneracy `P1 = P2`, which
  `agentC_work/CERT_second_door.py` proves unreachable.
* Anyone re-deriving hardness from "the curve is not secp256k1" is building on the a2 = 0 error.
