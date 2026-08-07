#!/usr/bin/env python3
"""Agent P: symbolic backward expansion to extract every 'stage' and its law."""
import pickle,sys,json
from collections import defaultdict,Counter
sys.set_int_max_str_digits(10**7)
sys.setrecursionlimit(100000)
W='/home/user/integer_solver/solve_lab/agentP_work/'
D=pickle.load(open(W+'model4.pkl','rb')); AP=D['AP']
S=pickle.load(open(W+'slp.pkl','rb')); topo=S['topo']; outof=S['outof']; definedat=S['definedat']
NV=38748
g=[0]*NV
for k,v in json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')).items(): g[int(k[2:])]=int(v)
P=115792089237316195423570985008687907853269984665640564039457584007908834671663
Q=97553848499418123410591666447050222001188385549510401465815187079080512838891
PV={x for x in range(NV) if g[x]==P}; QV={x for x in range(NV) if g[x]==Q}
atom_of={}                       # var -> atom index that defines it
for x,i in definedat.items(): atom_of[x]=topo[i]

# symbolic: polynomial over "leaf symbols" = free vars (+ special P,Q markers)
CAP=400
def pmul(a,b):
    r=defaultdict(int)
    for m1,c1 in a.items():
        for m2,c2 in b.items(): r[tuple(sorted(m1+m2))]+=c1*c2
    return {k:v for k,v in r.items() if v}
def padd(a,b):
    r=dict(a)
    for m,c in b.items():
        r[m]=r.get(m,0)+c
        if not r[m]: del r[m]
    return r

memo={}
def expand(x, depth=0):
    if x in memo: return memo[x]
    if x in PV: return {('P',):1}
    if x in QV: return {('Q',):1}
    if x not in atom_of or depth>40: return {(x,):1}
    ap=AP[atom_of[x]]
    # solve co*x + rest = 0  -> x = -rest/co
    co=0; rest={}
    for m,c in ap.items():
        if m==(x,): co+=c; continue
        t={():c}
        for y in m: t=pmul(t,expand(y,depth+1))
        rest=padd(rest,t)
    if co==0: return {(x,):1}
    out={}
    for m,c in rest.items():
        if c % co: return {(x,):1}     # non-integral -> treat as opaque
        out[m]=-c//co
    if len(out)>CAP: out={(x,):1}
    memo[x]=out
    return out

def sig(poly):
    """shape signature: relabel symbols by first occurrence, keep coeffs"""
    ren={}; terms=[]
    for m,c in sorted(poly.items(),key=lambda z:(-len(z[0]),z[0])):
        mm=[]
        for s in m:
            if s not in ren: ren[s]='s%d'%len(ren) if not isinstance(s,str) else s
            mm.append(ren[s])
        terms.append((tuple(sorted(mm)),c))
    return tuple(sorted(terms)), ren

if __name__=='__main__':
    # find Q-gates: atoms containing a Q-alias var and >=1 other var
    qg=[]
    for i,a in enumerate(topo):
        vs=set()
        for m in AP[a]: vs.update(m)
        if vs & QV and len(vs)>1: qg.append(i)
    print("atoms with Q-alias and other vars:",len(qg))
    shapes=Counter()
    for i in qg:
        a=topo[i]
        shapes[tuple(sorted((tuple(sorted(('Q' if x in QV else 'v') for x in m)),c if abs(c)<10**12 else 'BIG') for m,c in AP[a].items()))]+=1
    for s,c in shapes.most_common(10): print(c,s)
