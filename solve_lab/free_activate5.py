#!/usr/bin/env python3
"""CLEAN escape: the 220-var wire is a QUIET free parameter (uniform wire=V breaks
0 extra atoms). So freeze the WHOLE wire uniformly=V and the slack vars, activate the
slack via partners x_30077/x_31807/x_31302 (now a1813/a1815/a1816 are SATISFIED since
x_38215=x_30077*V etc. match), hold the twist, and count violations.
At x_12779=2 (x_14402=-1), V=1, knob t:
  x_24026 = 321447*t ; x_30077 = -t  (=> a1813 holds)
  x_9770  = x_35186 + 642894*t  (9770 gap absorbed by a1817 slack)
  x_27116 = (x_17728-x_1642)/2 ; x_31807 = -x_27116  (=> a1815 holds, x_3183=x_17728)
  x_31302 = 6033033*(x_18274 - x_9770)  (=> a1817/a1816 hold, twist terms cancel in Q40782)
Scan t for minimum violations."""
import json, time
from confluent_eval5 import build5, make_forward
from slack_active import viol_atoms
from propagate import atom_vars, NVARS

def main():
    t0=time.time()
    A,kind,info,seq0,bestval,ncyc=build5()
    order=json.load(open('eval_order.json'))['order']
    defset=set(v for v in kind if kind[v]!='const')
    seq=[v for v in order if v in defset and v not in (9770,3183)]
    seq+=[v for v in (9770,3183) if v in defset]
    seq+=[v for v in defset if v not in set(order) and v not in (9770,3183)]
    solve=make_forward(kind,info,seq,bestval)
    control=json.load(open('control_bits.json'))
    # wire class
    parent={}
    def f(x):
        parent.setdefault(x,x)
        while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
        return x
    def u(a,b):
        ra,rb=f(a),f(b)
        if ra!=rb: parent[ra]=rb
    for poly in A:
        if len(poly)==2:
            (m1,c1),(m2,c2)=list(poly.items())
            if len(m1)==1 and len(m2)==1 and abs(c1)==abs(c2): u(m1[0],m2[0])
    wire=set(x for x in list(parent) if f(x)==f(15))
    # solver that freezes wire + slack vars, recomputes rest
    FREEZE=wire | {24026,27116,26977}
    seq2=[v for v in seq if v not in FREEZE]
    def run(val, frozen):
        for v in FREEZE: val[v]=frozen.get(v,val[v])
        for v in seq2:
            k=kind[v]
            if k=='gate':
                coef,terms=info[v]; rs=0
                for c,m in terms:
                    t=c
                    for x in m: t*=val[x]
                    rs+=t
                if coef and (-rs)%coef==0: val[v]=(-rs)//coef
            elif k=='load':
                bit,cbx,lt=info[v]
                if val[bit]==0: val[v]=0
                else:
                    rest=0
                    for c,m in lt:
                        t=c
                        for x in m: t*=(1 if x==bit else val[x])
                        rest+=t
                    num=-rest; den=cbx*val[bit]
                    if den and num%den==0: val[v]=num//den
            elif k=='div':
                c,uu,rest=info[v]; rs=0
                for cc,m in rest:
                    t=cc
                    for x in m: t*=val[x]
                    rs+=t
                den=c*val[uu]
                if den and (-rs)%den==0: val[v]=(-rs)//den
                elif den==0: val[v]=0
        return val

    st=5
    def rnd():
        nonlocal st; st=(st*6364136223846793005+1442695040888963407)&((1<<64)-1); return st>>33
    # find x_12779=2 state with (x_17728-x_1642) even
    S=None
    for _ in range(1500):
        k=8+rnd()%30; cand=sorted(set(control[rnd()%len(control)] for _ in range(k)))
        v=solve(list(bestval),cand)
        if v[12779]==2 and v[18274]!=0 and (v[17728]-v[1642])%2==0: S=cand; v1=v; break
    if S is None: print("no state"); return
    print(f"state k={len(S)}; x_18274!=0, x_17728-x_1642 even ({time.time()-t0:.0f}s)", flush=True)
    V=1; x27116=(v1[17728]-v1[1642])//2
    best=None
    for t in [0,1,-1,2,-2,10,-10,1000,-1000, (v1[18274]-v1[35186])//642894]:
        x24026=321447*t; x9770=v1[35186]+642894*t
        x30077=-t; x31807=-x27116//V if x27116%V==0 else -x27116; x31302=6033033*(v1[18274]-x9770)//V
        val=list(v1)
        for w in wire: val[w]=V
        frozen={24026:x24026,27116:x27116,26977:V*x31302}
        val[30077]=x30077; val[31807]=x31807; val[31302]=x31302
        v2=run(val,frozen)
        tw=(v2[3183]==v2[17728]); a1817=(v2[26977]==6033033*(v2[18274]-v2[9770]))
        bad=viol_atoms(A,v2)
        if best is None or len(bad)<best[0]: best=(len(bad),t,sorted(bad))
        print(f"  t={t}: x3183==x17728 {tw}, a1817 {a1817}, x24026!=0 {v2[24026]!=0}, violated {len(bad)}: {sorted(bad)[:14]}", flush=True)
    print(f"\nbest: {best[0]} violated at t={best[1]}: {best[2][:16]} ({time.time()-t0:.0f}s)", flush=True)
    if best[0]==0: print("*** SOLVED ***")

if __name__=='__main__':
    main()
