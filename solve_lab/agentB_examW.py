#!/usr/bin/env python3
"""Examine wire1correct.json: failing set (exact Z), wire values, and for the 13 unpacking eqs
whether they depend on handles (free inputs) or only on wire members. Check L2 mod 6672769."""
import json, pickle
from collections import defaultdict
from agentB_setup import load, Env, p, NVARS

data = load(); env = Env(data)
wire = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/agentB_wire.pkl','rb'))['wire']
wireset = set(wire); freeset = env.freeset
sol = {int(k[2:]): int(v) for k, v in json.load(open('wire1correct.json')).items()}
valz = [0]*NVARS
for k, v in sol.items(): valz[k] = v
def root_int(i):
    s = 0
    for m, c in env.root_poly[i].items():
        t = c
        for v in m: t *= valz[v]
        s += t
    return s
fails = [i for i in range(len(env.root_poly)) if root_int(i) != 0]
UNPACK = set([8429,11166,11915,12594,23869,25313,26785,31400,32300,36106,36767,37257,37666])
CORE = set([2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892])
print(f"[W] wire1correct failing={len(fails)}: core={sorted(set(fails)&CORE)}")
print(f"    unpack={sorted(set(fails)&UNPACK)}  other={sorted(set(fails)-CORE-UNPACK)}")

# wire values
wv = set(valz[v] for v in wire)
print(f"[W] wire distinct values: {len(wv)}; x_26064={valz[26064]} (=1?{valz[26064]==1}) x_5101={valz[5101]}")

# quotient handles
for h in [30317, 2936, 5146]:
    print(f"[W] handle x_{h}={valz[h]}")

# does each unpacking eq depend on any FREE input (handle)?  Use eqvars restricted to freeset.
eqvars = env.eqvars
for i in sorted(UNPACK):
    fv = eqvars[i] & freeset
    wv2 = eqvars[i] & wireset
    print(f"[W] unpack eq {i}: free-inputs in eq={len(fv)} {sorted(fv)[:6]}, wire in eq={len(wv2)}, resid={root_int(i)%p if root_int(i)%p<10**6 else root_int(i)%p}")

# L2 = 10159099*S + 6926539*T ; S=x_35389, T=x_6671 ; L2=x_25739
print(f"[W] L2=x_25739={valz[25739]}  mod 6672769={valz[25739]%6672769}")
print(f"[W] S=x_35389 mod p={valz[35389]%p if valz[35389]%p<10**7 else '...'}, T=x_6671 mod p={valz[6671]%p if valz[6671]%p<10**7 else '...'}")
print(f"[W] x_3558={valz[3558]}, x_24908={valz[24908]}, x_16742={valz[16742]}")
