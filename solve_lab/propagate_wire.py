#!/usr/bin/env python3
"""Test the '256 freedom' hypothesis: does (control bits + wire value V) determine
the ~7053 free partner vars via propagation through the atoms? Merge the 220-var wire
-> x_15 and the twist x_17728->x_3183. Seed a bit-setting + x_15=V + pins. Propagate
over GF(P) (fast, exact for determination): atom with exactly ONE unknown (linear or
product/div solvable) -> assign. Report how many of the 38748 vars get determined and
contradictions. If nearly all determine, freedom ~= bits+V (searchable); if it stalls
low, the escape is genuinely high-dimensional."""
import json, time, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
from modp import P, inv

def main():
    t0=time.time()
    A0=load_atoms()
    control=list(json.load(open('control_bits.json')))
    # wire union-find (identity/negation atoms)
    parent={}
    def f(x):
        parent.setdefault(x,x)
        while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
        return x
    def u(a,b):
        ra,rb=f(a),f(b)
        if ra!=rb: parent[ra]=rb
    for poly in A0:
        if len(poly)==2:
            (m1,c1),(m2,c2)=list(poly.items())
            if len(m1)==1 and len(m2)==1 and abs(c1)==abs(c2): u(m1[0],m2[0])
    r15=f(15); wire=set(x for x in list(parent) if f(x)==r15)
    MERGE={x:15 for x in wire}; MERGE[17728]=3183   # twist merge
    def rep(v): return MERGE.get(v,v)

    # remap atoms mod P
    atoms=[]
    for poly in A0:
        out=defaultdict(int)
        for m,c in poly.items():
            nm=tuple(sorted(rep(x) for x in m)); out[nm]=(out[nm]+c)%P
        atoms.append({m:c for m,c in out.items() if c})

    # seed: bit-setting + x_15=V + pins-from-atoms
    st=int(sys.argv[1]) if len(sys.argv)>1 else 999
    V=int(sys.argv[2]) if len(sys.argv)>2 else 1
    def rnd():
        nonlocal st; st=(st*6364136223846793005+1442695040888963407)&((1<<64)-1); return st>>33
    val=[None]*NVARS
    nb=0
    for b in control:
        val[rep(b)]=rnd()&1; nb+=1
    bits=set()
    val[15]=V%P
    var_atoms=defaultdict(list)
    for a,poly in enumerate(atoms):
        for v in set().union(*[set(m) for m in poly]) if poly else []: var_atoms[v].append(a)
    contra=[0]
    def assign(v,x):
        x%=P
        if val[v] is not None:
            if val[v]!=x: contra[0]+=1
            return
        val[v]=x
    wl=deque(range(len(atoms)))
    def subst(poly):
        out=defaultdict(int)
        for m,c in poly.items():
            cc=c; nm=[]
            for v in m:
                if val[v] is not None: cc=(cc*val[v])%P
                else: nm.append(v)
            out[tuple(sorted(nm))]=(out[tuple(sorted(nm))]+cc)%P
        return {m:c for m,c in out.items() if c}
    passes=0
    while wl:
        passes+=1
        if passes>60: break
        nxt=deque()
        changed=False
        for a in range(len(atoms)):
            r=subst(atoms[a])
            uv=set().union(*[set(m) for m in r]) if r else set()
            if len(uv)==0:
                if r.get((),0)%P: contra[0]+=1
                continue
            if len(uv)==1:
                v=next(iter(uv)); c0=c1=c2=0
                for m,c in r.items():
                    if len(m)==0: c0=(c0+c)%P
                    elif len(m)==1: c1=(c1+c)%P
                    else: c2=(c2+c)%P
                if c2==0 and c1: assign(v,(-c0*inv(c1))%P); changed=True
                elif c2==0 and c1==0:
                    if c0%P: contra[0]+=1
        if not changed: break
    ndet=sum(1 for x in val if x is not None)
    ndetfree=sum(1 for v in range(NVARS) if val[v] is not None and rep(v)==v)
    print(f"wire size {len(wire)}; seed bits={len(bits)}, V={V}", flush=True)
    print(f"propagation: {passes} passes, determined {ndet}/{NVARS} rep-vars, contradictions {contra[0]} ({time.time()-t0:.0f}s)", flush=True)
    print(f"  (freedom test: if determined ~= NVARS-256, the escape is a (bits,V) search)", flush=True)

if __name__=='__main__':
    main()
