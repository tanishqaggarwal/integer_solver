#!/usr/bin/env python3
"""Diagnose a config: for each failing eq, flatten into atoms and print nonzero atoms + values.
Also print values of routing handles and control vars."""
import sys, json, ast
import agentD_harness as H
C1,C2=H.C1,H.C2
mode=sys.argv[1] if len(sys.argv)>1 else '10'
CONST={30213:C2, 22162:C1, 24468:C1, 18956:C2}
a7=24601; a34=2081
if mode=='11': ov={a7:1,a34:1,**CONST}
elif mode=='10': ov={a7:1,**CONST}
elif mode=='01': ov={a34:1,**CONST}
elif mode=='00': ov={**CONST}
elif mode=='00noC': ov={}
r=H.run_config(ov, want_val=True)
H.val=r['val']; H.ns['v']=H.val
print(f"mode {mode}: sat={r['satisfied']} nfail={r['nfail']} x15298={r['x_15298']} S0={r['S_is0']} T0={r['T_is0']}")
print("F:", r['F'])
# handle values
handles=[13682,37892,34243,32237,30213,22162,24468,18956,15298,11150,25739,37758,35605,4007,29804,
         30317,2936,5146,35389,6671,3558,29322,14853,12186,16742,24908,1326]
print("\n=== handle/control values (mod p shown if huge) ===")
for h in handles:
    v=H.val[h]
    tag='=C1' if v==C1 else '=C2' if v==C2 else '=-C1' if v==-C1 else '=-C2' if v==-C2 else ''
    vs=str(v)
    if len(vs)>25: vs=f"{vs[:12]}...({len(vs)}d) mod p={v%H.p}"
    print(f"  x_{h}: {vs} {tag}")
# per failing eq: flatten and show nonzero atoms
print("\n=== failing-eq atom breakdown ===")
for idx in r['F']:
    root=H.rootast(idx)
    terms=H.flat(root)
    nz=[]
    for t in terms:
        val=H.evn(t)
        if val!=0:
            src=ast.unparse(t)
            if len(src)>70: src=src[:70]+'...'
            nz.append((val,src))
    print(f"eq{idx}: {len(nz)} nonzero atoms of {len(terms)}")
    for val,src in nz[:6]:
        vs=str(val)
        if len(vs)>20: vs=f"({len(vs)}d, modp={val%H.p})"
        print(f"    {vs}: {src}")
