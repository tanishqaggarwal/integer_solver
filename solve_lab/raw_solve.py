#!/usr/bin/env python3
"""Solve the twist by requiring the raw EQUATIONS to vanish (not each leaf atom).

Key insight (as in the prior instance): the equations constrain linear *combinations* of
the shared atoms; a perfect-square equation E^2=0 only needs its root E=0. So atoms need
not be individually zero — they may cancel. Holding the V0-pinned wire and all non-twist
vars fixed, every twist atom is linear in the free twist vars, so each failing equation
(square -> its root) becomes a linear constraint. Solve that linear system over Z (SNF)."""
import json, ast, re, sys
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS

base = {int(k[2:]): x for k, x in json.load(open('rebuilt_partial.json')).items()}
A = load_atoms()
# wire (220-class) held at V0
parent = {}
def f(x):
    parent.setdefault(x, x)
    while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
    return x
def u(a, b):
    ra, rb = f(a), f(b)
    if ra != rb: parent[ra] = rb
for poly in A:
    if len(poly) == 2:
        (m1,c1),(m2,c2)=list(poly.items())
        if len(m1)==1 and len(m2)==1 and abs(c1)==abs(c2): u(m1[0],m2[0])
classes=defaultdict(list)
for x in list(parent): classes[f(x)].append(x)
wire=set(max(classes.values(),key=len))

# unknown twist vars = every ADJUSTABLE var (small base, non-wire) in the failing equations.
# Everything else held at base; wire held at V0. Requiring EQUATIONS (not atoms) = 0.
FAIL=[56,133,2071,2683,4386,7254,8073,11009,13660,15299,16622,17726,19066,19656,19712,
      20452,22093,25480,28090,28653,31061,31138,32894,34517,34892,35089,35299,38629]
_lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
UNKS=set()
for i in FAIL:
    for z in set(int(m) for m in re.findall(r'x_(\d+)', _lines[i])):
        if z not in wire and abs(base.get(z,0))<10**6: UNKS.add(z)
UNK=sorted(UNKS); idx={v:i for i,v in enumerate(UNK)}; n=len(UNK)

def snf_solve(Am,b):
    m=len(Am); nn=len(Am[0]) if Am else 0
    M=[r[:] for r in Am]; bb=b[:]
    V=[[1 if i==j else 0 for j in range(nn)] for i in range(nn)]
    for piv in range(min(m,nn)):
        while True:
            best=None;found=None
            for i in range(piv,m):
                for j in range(piv,nn):
                    if M[i][j]!=0 and (best is None or abs(M[i][j])<best): best=abs(M[i][j]);found=(i,j)
            if found is None: break
            i0,j0=found
            if i0!=piv: M[piv],M[i0]=M[i0],M[piv]; bb[piv],bb[i0]=bb[i0],bb[piv]
            if j0!=piv:
                for r in M: r[piv],r[j0]=r[j0],r[piv]
                for r in V: r[piv],r[j0]=r[j0],r[piv]
            p=M[piv][piv]; done=True
            for i in range(m):
                if i!=piv and M[i][piv]!=0:
                    q=M[i][piv]//p
                    for j in range(nn): M[i][j]-=q*M[piv][j]
                    bb[i]-=q*bb[piv]
                    if M[i][piv]!=0: done=False
            for j in range(nn):
                if j!=piv and M[piv][j]!=0:
                    q=M[piv][j]//p
                    for i in range(m): M[i][j]-=q*M[i][piv]
                    for i in range(nn): V[i][j]-=q*V[i][piv]
                    if M[piv][j]!=0: done=False
            if done: break
    y=[0]*nn
    for i in range(min(m,nn)):
        d=M[i][i]
        if d==0: continue
        if bb[i]%d!=0: return None
        y[i]=bb[i]//d
    for i in range(m):
        if all(z==0 for z in M[i]) and bb[i]!=0: return None
    return [sum(V[i][j]*y[j] for j in range(nn)) for i in range(nn)]

# substitute non-UNK vars (wire->V0, others->base) at the AST level, expand to poly in UNK
V0=base[next(iter(wire))]
def expand(node):
    if isinstance(node,ast.Constant): return {(): node.value}
    if isinstance(node,ast.Name):
        vi=int(node.id[2:])
        if vi in UNKS: return {(vi,):1}
        return {(): (V0 if vi in wire else base.get(vi,0))}
    if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):
        return {m:-c for m,c in expand(node.operand).items()}
    if isinstance(node,ast.BinOp):
        a=expand(node.left); b=expand(node.right)
        if isinstance(node.op,ast.Add):
            r=dict(a)
            for m,c in b.items(): r[m]=r.get(m,0)+c
            return {m:c for m,c in r.items() if c}
        if isinstance(node.op,ast.Sub):
            r=dict(a)
            for m,c in b.items(): r[m]=r.get(m,0)-c
            return {m:c for m,c in r.items() if c}
        if isinstance(node.op,ast.Mult):
            r={}
            for m1,c1 in a.items():
                for m2,c2 in b.items():
                    m=tuple(sorted(m1+m2)); r[m]=r.get(m,0)+c1*c2
            return {m:c for m,c in r.items() if c}
    raise ValueError

def strip_sq(node):
    while isinstance(node,ast.BinOp) and isinstance(node.op,ast.Mult):
        a,b=node.left,node.right
        ca=isinstance(a,ast.Constant) or (isinstance(a,ast.UnaryOp) and isinstance(a.operand,ast.Constant))
        cb=isinstance(b,ast.Constant) or (isinstance(b,ast.UnaryOp) and isinstance(b.operand,ast.Constant))
        if ca and not cb: node=b
        elif cb and not ca: node=a
        elif (not ca) and (not cb):
            if ast.unparse(a)==ast.unparse(b): node=a
            else: break
        else: break
    return node

lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
rows=[]; bs=[]; nonlin=0; used=0
for L in lines:
    ids=set(int(m) for m in re.findall(r'x_(\d+)',L))
    if not (ids & UNKS): continue
    used+=1
    lhs=L.rsplit('=',1)[0]
    node=strip_sq(ast.parse(lhs,mode='eval').body)
    poly=expand(node)
    coef=[0]*n; const=0; ok=True
    for m,c in poly.items():
        us=[z for z in m if z in idx]
        if len(us)==0: const+=c
        elif len(us)==1 and len(m)==1: coef[idx[us[0]]]+=c
        else: ok=False; break   # UNK*UNK or UNK*held-nonzero handled below
    if not ok:
        # linearize: hold all but first UNK factor at base
        coef=[0]*n; const=0
        for m,c in poly.items():
            us=[z for z in m if z in idx]
            if not us:
                const+=c
            else:
                keep=us[0]; t=c
                for z in m:
                    if z==keep: continue
                    t*=(V0 if z in wire else base.get(z,0)) if z not in idx else base.get(z,0)
                coef[idx[keep]]+=t
        nonlin+=1
    rows.append(coef); bs.append(-const)
print(f"{used} equations touch twist vars; {len(rows)} constraints ({nonlin} linearized), {n} unknowns", flush=True)
x=snf_solve(rows,bs)
if x is None:
    print("NO integer solution"); sys.exit(0)
cand=dict(base)
for v in UNK: cand[v]=x[idx[v]]
json.dump({f"x_{i}":cand.get(i,0) for i in range(NVARS)}, open('raw_solved.json','w'))
print("wrote raw_solved.json", flush=True)
