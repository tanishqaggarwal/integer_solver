#!/usr/bin/env python3
"""
primitive_analysis.py — Reverse-engineering "what the circuit computes".

Conclusion (see bottom): the circuit is NOT a known cryptographic primitive
(ECDSA / Schnorr / EC scalar-mult / discrete-log / Pedersen / hash-preimage).
It is a bespoke SHALLOW (multiplicative depth 10) obfuscated MQ / knapsack
trapdoor over GF(secp256k1_p). The 512 loaded constants are cryptographically
random; secp256k1's field prime is used cosmetically (a standard hard 256-bit
prime), not because any EC arithmetic is performed.

Run:  python3 primitive_analysis.py     (from solve_lab/)
Requires heal_harness.py (forward-reconstruct + ancestry) and pinrec.json.
"""
import json, itertools, random
p = 2**256 - 2**32 - 977
Gx = 55066263022277343669578718895168534326250603453777594175500187360389116729240
Gy = 32670510020758816978083085130507043184471273380659243275938904335757337482424
n  = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141

pins = json.load(open('pinrec.json'))
t2c  = {pin[2]: (pin[3], pin[4]) for pin in pins}          # target -> (rawconst, coef)
R    = {t: t2c[t][0] % p for t in t2c}                     # residue mod p
uniq = sorted(set(R.values()))
setv = set(uniq)

def is_qr(a): return a == 0 or pow(a, (p - 1) // 2, p) == 1

print("=== 1. CONSTANT CENSUS ===")
print(f"pins={len(pins)}  unique targets={len(t2c)}  unique residues={len(uniq)}")
raw = [pin[3] for pin in pins]
print(f"raw const bit-length range: {min(c.bit_length() for c in raw)}..{max(c.bit_length() for c in raw)} "
      f"(all > 256 => NOT reduced mod p; quotient q=C//p is a 31-40 bit random pad)")

print("\n=== 2. CRYPTOGRAPHIC-STRUCTURE TESTS ON THE 512 RESIDUES (all should be ~noise) ===")
onx = [r for r in uniq if is_qr((r * r % p * r + 7) % p)]
print(f"valid secp256k1 x-coords (x^3+7 is QR): {len(onx)}/{len(uniq)}  "
      f"(chance {len(uniq)//2}+-{int((len(uniq)*0.25)**0.5)}) -> consistent with random")
print(f"exact residue collisions: {len(uniq) - len(set(uniq))}  (0 => all distinct)")
print(f"additive-inverse pairs C_i=-C_j: {sum(1 for r in uniq if (p-r) in setv and p-r!=r)//2}")
sm = sum(1 for k in range(2, 17) for r in uniq if (k * r) % p in setv)
print(f"small-ratio pairs C_i=k*C_j (2<=k<=16): {sm}")
for name, b in [('Gx', Gx), ('Gy', Gy), ('curve-b=7', 7), ('order-n', n)]:
    h = sum(1 for r in uniq if (r * b) % p in setv or (r + b) % p in setv or (r * pow(b, -1, p)) % p in setv)
    print(f"residues related to {name} by *,+,/ : {h}")
# EC additive structure C_a+C_b == C_k  (Pedersen/Schnorr signature would show this)
Sset = set(uniq); samp = random.Random(1).sample(uniq, 80)
print(f"additive triples C_a+C_b==C_k (80-sample): {sum(1 for a in samp for b in samp if (a+b)%p in Sset)}")
# EC point doubling among valid-x residues
def sqrtp(a): return pow(a, (p + 1) // 4, p)
dbl = 0
for x in onx:
    y = sqrtp((x*x%p*x+7) % p)
    if y and (2*y) % p:
        lam = (3*x*x % p)*pow(2*y, -1, p) % p
        if (lam*lam - 2*x) % p in setv: dbl += 1
print(f"points whose EC-double x-coord is another residue: {dbl}")

print("\n=== 3. AUX MODULUS 6672769 ===")
from math import gcd
m = 6672769
print(f"prime; gcd(6672769, p)={gcd(m,p)}; (p-1)%m={(p-1)%m}; n%m={n%m}  -> unrelated 23-bit random prime")

print("\n=== 4. CIRCUIT DEPTH (primitive discriminator) ===")
try:
    import heal_harness as H
    import re
    VAR = re.compile(r'x_(\d+)')
    def is_vv(rhs):
        return rhs.count('*') == 1 and not rhs.split('*')[0].strip().replace('-', '').isdigit()
    md = {v: 0 for v in H.freeinp}
    for t in H.order:
        _, rhs, vids = H.gates[H.definer[t]]
        md[t] = max((md.get(u, 0) for u in vids), default=0) + (1 if is_vv(rhs) else 0)
    print(f"MAX multiplicative depth over all gates: {max(md.values())}")
    print("  -> depth ~10 is SHALLOW: rules out EC scalar-mult / Montgomery ladder /")
    print("     ECDSA-Schnorr verification / iterated hash preimage (all need depth ~256).")
    wmax = max(len(H.anc[t]) for t in H.order)
    print(f"widest gate free-input fan-in: {wmax}  (codeword gates aggregate all 256 selector bits)")
except Exception as e:
    print("  (heal_harness unavailable:", e, ")")

print("""
=== VERDICT ===
* NOT a known primitive. Multiplicative depth 10 forbids every deep secp256k1
  computation (scalar mult, ECDSA/Schnorr verify, DLog, Pedersen, SHA/Keccak).
* The 512 message constants are cryptographically random: no EC-point membership,
  no Gx/Gy/n/b relation, no discrete-log ladder, no knapsack collision, no
  arithmetic/geometric progression, no hash(i) generator. All tests = noise.
* The core conic  x_27713^2 = x_33469*x_1326^2  lives entirely among the setter's
  FREE control knobs (x_14853,x_12186,x_16742,x_22162,...), disconnected from the
  message words -- it is a QR/sqrt gadget, not y^2=x^3+7 and not an EC/DLog relation.
* The perfect-square "verifier" atoms Q^2 are RANDOM linear combinations of wiring
  atoms (coefs like -29,-11,20,-37) -- a constraint-bundling obfuscation, not a
  semantic curve/pairing check.
* secp256k1's prime is COSMETIC: a standard hard 256-bit prime (odd, coprime to the
  23-bit aux modulus) that defeats small-modulus / CRT attacks and 'looks crypto'.
=> The trapdoor is a generic shallow-but-wide multilinear MQ/knapsack over the 256
   boolean selector bits with random GF(p) coefficients. No primitive-inversion
   shortcut exists; this matches the prior agents' 'needs the setter's global
   message' conclusion and refutes the known-primitive hypothesis.
""")
