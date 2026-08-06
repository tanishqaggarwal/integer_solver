#!/usr/bin/env python3
"""Diagnose an explicit override (passed as JSON on argv[1]).
e.g. python3 agentD_diag2.py '{"24601":1,"24468":"C1","18956":"C2"}'  (C1/C2 as strings)."""
import sys, json, ast
import agentD_harness as H
C1,C2=H.C1,H.C2
raw=json.loads(sys.argv[1])
ov={}
for k,v in raw.items():
    kk=int(k)
    if v=='C1': ov[kk]=C1
    elif v=='C2': ov[kk]=C2
    else: ov[kk]=int(v)
r=H.run_config(ov, want_val=True)
H.val=r['val']; H.ns['v']=H.val
print(f"override keys: {sorted(ov)}")
print(f"sat={r['satisfied']} nfail={r['nfail']} x15298={r['x_15298']} S0={r['S_is0']} T0={r['T_is0']} core={r['core_fail']} noncore={r['noncore_fail']}")
print("F:", r['F'])
print("S mod p =", r['S_modp'])
print("T mod p =", r['T_modp'])
# atom breakdown for each failing eq
print("\n=== failing-eq atom breakdown ===")
atom_count={}
for idx in r['F']:
    root=H.rootast(idx)
    terms=H.flat(root)
    nz=[]
    for t in terms:
        v=H.evn(t)
        if v!=0:
            src=ast.unparse(t)
            # strip leading integer coeff for grouping
            key=src
            nz.append((v,src))
            atom_count[src]=atom_count.get(src,0)+1
    incore = idx in H.CORESET
    print(f"eq{idx}{'[core]' if incore else '[NONcore]'}: {len(nz)} nz atoms")
    for v,src in nz[:6]:
        vs=str(v)
        if len(vs)>18: vs=f"({len(vs)}d,modp={v%H.p})"
        print(f"    {vs}: {src[:80]}")
print("\n=== atom frequency (which gate-atoms recur) ===")
for src,c in sorted(atom_count.items(), key=lambda x:-x[1])[:15]:
    print(f"  x{c}: {src[:90]}")
