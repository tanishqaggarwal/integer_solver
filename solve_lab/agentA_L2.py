#!/usr/bin/env python3
"""M2 subtask: L2 = x_25739 = 10159099*S+6926539*T (exact). Analyze L2 mod 6672769 as a
function of x_3558 (message-controlled scalar) and message bits; solve L2 ≡ 0 mod 6672769."""
import json, os
from agentA_harness import (p, load_solution, forward, backward_cone, freeinp, eqcode,
                            eqvars, NEQ)
M = 6672769   # prime
base = load_solution('best/new_instance_partial_39013.json'); forward(base)
boolset = set(json.load(open('boolbits.json'))['boolvars'])
allbits = sorted((backward_cone(35389)[1] | backward_cone(6671)[1]) & boolset)

def L2(v): return v[25739]
def count(v):
    ns = {'__builtins__': {}, 'v': v}
    return set(i for i in range(NEQ) if eval(eqcode[i], ns) != 0)

# --- L2 mod M as quadratic in y = x_3558 (verify) ---
# S = x_33469*x_29322^2 - x_3558^2 ; T = x_27713*x_29322 - x_3558*x_1326 ; x_3558 = x_24908 - x_16742
a = base[33469]; d = base[29322]; u = base[27713]; w = base[1326]; k = base[16742]; X = base[24908]
def L2_formula(y):
    S = a*d*d - y*y
    T = u*d - y*w
    return 10159099*S + 6926539*T
y0 = base[3558]
print(f"L2 formula check: {L2_formula(y0) == base[25739]}  (y0=x_3558)")
# quadratic coeffs mod M: L2 = C0 - 10159099*y^2 - 6926539*w*y
C0 = (10159099*a*d*d + 6926539*u*d) % M
c2 = (-10159099) % M
c1 = (-6926539*w) % M
print(f"L2 mod {M} = {C0} + {c1}*y + {c2}*y^2   (y=x_3558 mod {M})")
print(f"  check at y0: {(C0 + c1*y0 + c2*y0*y0) % M} == {base[25739]%M} ? {(C0+c1*y0+c2*y0*y0)%M==base[25739]%M}")

# --- solve quadratic c2*y^2 + c1*y + C0 ≡ 0 mod M ---
def sqrt_mod(n, q):
    n %= q
    if n == 0: return [0]
    if pow(n, (q-1)//2, q) != 1: return []
    if q % 4 == 3:
        r = pow(n, (q+1)//4, q); return sorted({r, q-r})
    # Tonelli-Shanks
    import random
    qq = q-1; s = 0
    while qq % 2 == 0: qq//=2; s+=1
    z = 2
    while pow(z, (q-1)//2, q) != q-1: z += 1
    mR = pow(n, (qq+1)//2, q); t = pow(n, qq, q); c = pow(z, qq, q); mm = s
    while t != 1:
        i = 0; tt = t
        while tt != 1: tt = tt*tt % q; i += 1
        b = pow(c, 1 << (mm-i-1), q); mR = mR*b % q; c = b*b % q; t = t*c % q; mm = i
    return sorted({mR, q-mR})
# c2 y^2 + c1 y + C0 = 0 -> y = (-c1 ± sqrt(c1^2-4 c2 C0)) / (2 c2)
disc = (c1*c1 - 4*c2*C0) % M
roots_disc = sqrt_mod(disc, M)
inv2c2 = pow((2*c2) % M, M-2, M)
yroots = sorted({((-c1 + rd) * inv2c2) % M for rd in roots_disc})
print(f"\ndiscriminant QR? {len(roots_disc)>0}; y = x_3558 mod {M} roots for L2≡0: {yroots}")
for yr in yroots:
    print(f"  y={yr}: L2 formula mod M = {L2_formula(yr) % M}  (want 0)")

# --- single-bit sensitivity of L2 mod M (message dependence) ---
print(f"\n=== L2 mod {M} single-bit sensitivities ===")
L20 = base[25739] % M
nz = 0
for b in allbits:
    v = base[:]; v[b] = 1 - base[b]; forward(v)
    if (v[25739] - base[25739]) % M != 0: nz += 1
print(f"bits changing L2 mod {M}: {nz}/{len(allbits)}  (baseline L2 mod M = {L20})")

# --- my 39021 config already has L2=0; confirm and copy as candidate ---
if os.path.exists('best_agentA_39021.json'):
    v = load_solution('best_agentA_39021.json'); forward(v)
    F = count(v)
    print(f"\nbest_agentA_39021: L2 mod {M} = {v[25739]%M}  x_25739 = {v[25739]}  wiring {NEQ-len(F)}/{NEQ}")
    if v[25739] % M == 0:
        json.dump({f"x_{i}": v[i] for i in range(len(v)) if v[i] != 0}, open('agentA_L2zero.json', 'w'))
        print("SAVED agentA_L2zero.json (L2 ≡ 0 mod 6672769, 39021 wiring)")
