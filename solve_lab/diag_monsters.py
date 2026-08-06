#!/usr/bin/env python3
"""Verify all 20 core eqs are integer combos of M1,M2,M3. Then decompose M1,M2,M3 into free inputs.
Critically: check whether those free inputs appear in the 39,013 currently-satisfied equations."""
import json, re, ast, sys
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
sys.setrecursionlimit(1000000)
p=2**256-2**32-977
sol={int(k[2:]):v for k,v in json.load(open('best/new_instance_partial_39013.json')).items()}
val=[0]*NVARS
for v,x in sol.items():
    if v<NVARS: val[v]=x
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
gate_defs={}
for t,rhs,vids in gates: gate_defs[t]=(rhs,vids)
gate_out=set(gate_defs); freeinp=set(v for v in range(NVARS) if v not in gate_out)
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
CORE=[2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892]
VAR=re.compile(r'x_(\d+)')
def evn(node):
    if isinstance(node,ast.Constant): return node.value
    if isinstance(node,ast.Name): return val[int(node.id[2:])]
    if isinstance(node,ast.UnaryOp): return -evn(node.operand)
    a=evn(node.left); b=evn(node.right)
    return a+b if isinstance(node.op,ast.Add) else a-b if isinstance(node.op,ast.Sub) else a*b
def inner(lhs):
    node=ast.parse(lhs,mode='eval').body
    while isinstance(node,ast.BinOp) and isinstance(node.op,ast.Mult):
        a,b=node.left,node.right
        ca=isinstance(a,ast.Constant) or (isinstance(a,ast.UnaryOp) and isinstance(a.operand,ast.Constant))
        cb=isinstance(b,ast.Constant) or (isinstance(b,ast.UnaryOp) and isinstance(b.operand,ast.Constant))
        if ca and not cb: node=b
        elif cb and not ca: node=a
        elif ast.unparse(a)==ast.unparse(b): node=a
        else: break
    return node
# define M1,M2,M3 by evaluating with x_15298=1 and current vals
M1 = val[15298]*val[11150] + val[4007]
M2 = val[15298]*val[25739] - 6672769*val[29804]
M3 = 537773*(val[15298]*val[37758]) - val[35605]
print(f"M1 = {M1}  (bits {M1.bit_length()})")
print(f"M2 = {M2}  (bits {M2.bit_length()})")
print(f"M3 = {M3}  (bits {M3.bit_length()})")
print(f"M1 mod p = {M1%p}")
print(f"M2 mod p = {M2%p}")
print(f"M3 mod p = {M3%p}")
print(f"M1//p={M1//p}, M2//p={M2//p}, M3//p={M3//p}")
# Solve for coefficients a,b,c s.t. E_i = a*M1+b*M2+c*M3 -- by evaluating E_i with perturbations
# Instead: E_i is exact integer; express via the three by GF: build integer least via matching.
# Simpler: recompute E_i, then greedily reduce by M3,M1,M2 (they have distinct bit-lengths / gcd)
print("\nPer-eq combos (E = a*M1 + b*M2 + c*M3):")
import itertools
for i in CORE:
    E=evn(inner(lines[i].rsplit('=',1)[0]))
    # solve small integer combo by trying a,b,c in range
    found=None
    for a in range(-40,41):
        r1=E-a*M1
        for c in range(-40,41):
            r2=r1-c*M3
            if M2!=0 and r2% M2==0:
                b=r2//M2
                if -40<=b<=40:
                    found=(a,b,c); break
        if found: break
    print(f"  eq {i}: {found}  {'OK' if found and found[0]*M1+found[1]*M2+found[2]*M3==E else 'MISMATCH E='+str(E)[:40]}")
