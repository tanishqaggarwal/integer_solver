#!/usr/bin/env python3
"""Task item 5: bit-flip search. From the (1,0) 39018 base, set each candidate boolvar=1
(one at a time) and re-run forward construction; record any that reduces fail count or drives
S/T residues toward 0. Candidates = boolvars in the free-input cones of the failing equations."""
import sys, json, time
import agentD_harness as H
C1,C2=H.C1,H.C2
p=H.p
boolvars=set(json.load(open('boolbits.json'))['boolvars'])
BASE={24601:1, 30213:0, 22162:0, 24468:C1, 18956:C2}  # (1,0) 39018
r0=H.run_config(BASE)
base_sat=r0['satisfied']; base_S=r0['S_modp']; base_T=r0['T_modp']
print(f"base (1,0) 39018: sat={base_sat} nfail={r0['nfail']} S%p!=0={base_S!=0} T%p!=0={base_T!=0}")
Ffail=r0['F']
# candidate boolvars: free-input ancestors of the failing equations, intersect boolvars, exclude base
cand=set()
for i in Ffail:
    for v in H.eqvars[i]:
        cand |= H.anc.get(v, {v} if v in H.freeinp else set())
cand &= boolvars
cand -= set(BASE)
cand = sorted(cand)
print(f"candidate boolvars in failing cones: {len(cand)}")
t0=time.time()
improved=[]; toward=[]
for k,b in enumerate(cand):
    r=H.run_config({**BASE, b:1})
    if r['satisfied']>base_sat:
        improved.append((r['satisfied'], b))
        print(f"  *** flip x_{b}=1: sat={r['satisfied']} (>{base_sat})!", flush=True)
        r2=H.run_config({**BASE, b:1}, want_val=True)
        json.dump({f"x_{i}":r2['val'][i] for i in range(H.NVARS)}, open(f"best_agentD_{r['satisfied']}.json",'w'))
    # track residue movement toward 0 (only meaningful if it changes S or T)
    if (r['S_modp']==0 and base_S!=0) or (r['T_modp']==0 and base_T!=0):
        toward.append((b, r['S_modp']==0, r['T_modp']==0, r['satisfied']))
    if k%40==0: print(f"  [{k}/{len(cand)}] x_{b}: sat={r['satisfied']} S0={r['S_is0']} T0={r['T_is0']} [{time.time()-t0:.0f}s]", flush=True)
print(f"\nIMPROVED flips (sat>{base_sat}): {improved}")
print(f"flips driving S or T to 0: {toward[:20]}")
# best-effort small combos among any flips that reduced nfail even if not above base
print(f"scanned {len(cand)} single flips in {time.time()-t0:.0f}s")
