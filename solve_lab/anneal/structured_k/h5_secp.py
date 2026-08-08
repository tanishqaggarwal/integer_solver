#!/usr/bin/env python3
"""h5_secp.py -- hypothesis 5: is this instance a re-coordinatised secp256k1 with
a *known* generator?  If  to_secp(G) == c*SECP_G  for a findable c, then
k = c^-1 * dlog(to_secp(T)) and the problem rebases onto the standard curve.

Subtlety: B/7 has SIX sixth roots mod p, so there are six F_p-isomorphisms onto
y^2 = x^3+7.  They differ by the order-6 automorphism (multiplication by a
primitive 6th root of unity in End(E)), so a small c under one root becomes
c*zeta^j under another.  We therefore test ALL SIX images of G and of T.
"""
import sys, os, subprocess, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from instance import p, n, B, G, T, PTS, add, mul, neg, sub, SECP_G
from sympy.ntheory.residue_ntheory import nthroot_mod

SK = os.path.join(HERE, 'sk')
BB = 26
GC = 1 << 22          # span 2^48 per target

# all sixth roots of B/7
target_val = B * pow(7, -1, p) % p
roots = sorted(set(nthroot_mod(target_val, 6, p, all_roots=True)))
print(f"B/7 has {len(roots)} sixth roots mod p")
assert all(pow(u, 6, p) == target_val for u in roots)

def iso(u, P):
    if P is None: return None
    return (P[0] * pow(u, -2, p) % p, P[1] * pow(u, -3, p) % p)

def on_secp(P):
    return (P[1] * P[1] - P[0] ** 3 - 7) % p == 0

targets = []
for j, u in enumerate(roots):
    Gi, Ti = iso(u, G), iso(u, T)
    assert on_secp(Gi) and on_secp(Ti), "isomorphism broken"
    targets.append((f'G_root{j}', Gi))
    targets.append((f'T_root{j}', Ti))

# the reverse direction: is SECP_G a small multiple of G on the instance curve?
targets_rev = [('SECPG_on_instance', None)]

# ---- planted control on the secp side ----
CONTROL_C = 987654321
def secp_mul(k, P):
    # y^2 = x^3 + 7 uses the same group law code (A == 0, B is not used)
    return mul(k, P)
ctrl_pt = secp_mul(CONTROL_C, SECP_G)
targets.append(('CONTROL', ctrl_pt))

if __name__ == '__main__':
    t0 = time.time()
    inp = "\n".join(f"{name} {P[0]:064x} {P[1]:064x}" for name, P in targets) + "\n"
    r = subprocess.run([SK, 'bsgsmulti', f"{SECP_G[0]:064x}", f"{SECP_G[1]:064x}",
                        str(BB), str(GC)], input=inp, capture_output=True, text=True)
    print(r.stderr.strip())
    M = 1 << BB
    results, found = [], []
    for line in r.stdout.split('\n'):
        f = line.split()
        if not f: continue
        if f[0] in ('MHIT', 'MNOHIT'):
            results.append({'target': f[1], 'hit': f[0] == 'MHIT'})
        elif f[0] == 'HITDATA':
            name = f[1]; i = int(f[2][2:]); jj = int(f[3][2:])
            for s in (1, -1):
                c = (i * M + s * jj) % n
                if secp_mul(c, SECP_G) == dict(targets)[name]:
                    found.append((name, c))
                    print(f"  {name} = {c} * SECP_G   (verified)")
    ctrl_ok = any(nm == 'CONTROL' and c == CONTROL_C for nm, c in found)
    print(f"control (c={CONTROL_C}): {'RECOVERED' if ctrl_ok else 'MISSED'}")

    span = M * GC
    print(f"\nsearched dlog base SECP_G over [0,{span}) = [0,2^{span.bit_length()-1}) "
          f"for {len(targets)} points")
    real = [(nm, c) for nm, c in found if nm != 'CONTROL']
    print("small dlogs found for instance points:", real if real else "NONE")

    # If we learned c with G_rootj = c*SECP_G we can rebase; check it solves.
    solved = None
    for nm, c in real:
        if nm.startswith('G_root'):
            j = int(nm[6:]); u = roots[j]
            # dlog of T under the same root: T_iso = k*G_iso = k*c*SECP_G
            print(f"  G maps to {c}*SECP_G under root {j}; k = dlog(T)/c")
        if nm.startswith('T_root'):
            j = int(nm[6:])
            print(f"  T maps to {c}*SECP_G under root {j}")

    json.dump({'babybits': BB, 'giantcount': GC, 'span': span,
               'span_bits': span.bit_length() - 1,
               'n_sixth_roots': len(roots), 'sixth_roots': [str(u) for u in roots],
               'targets': [t[0] for t in targets], 'results': results,
               'control_ok': ctrl_ok, 'found': [[nm, str(c)] for nm, c in found],
               'wall_seconds': round(time.time() - t0, 1)},
              open(os.path.join(HERE, 'h5_secp.json'), 'w'), indent=1)
    print(f"wall {time.time()-t0:.1f}s")
