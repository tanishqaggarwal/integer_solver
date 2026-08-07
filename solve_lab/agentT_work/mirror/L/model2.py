"""Complete tree model with NUMERIC pin extraction and NUMERIC per-node law verification."""
import sys, os, json, collections, pickle, re, time
F='/home/user/integer_solver/solve_lab/agentT_work/mirror/F'; sys.path.insert(0,F)
from fwd import Engine, NV
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
K=97553848499418123410591666447050222001188385549510401465815187079080512838891
defrhs={c[1]:c[2] for c in (E.cls[a] for a in E.order)}
H=pickle.load(open('handles.pkl','rb')); appearP=H['appearP']
D=pickle.load(open('ortree2.pkl','rb')); tree=D['tree']
NODE=pickle.load(open('nodes.pkl','rb')); OUT=pickle.load(open('outwires.pkl','rb'))
ridx=E.residx

def run(assign):
    v=[0]*NV
    for k,x in assign.items(): v[k]=x
    r=E.run(v)
    return v,r

leafnode={}
for n,N in NODE.items():
    for side,ch in (('va',N['a']),('vb',N['b'])):
        if tree[ch] is None: leafnode[ch]=(n,side)
live=[v for v in tree if tree[v] is None and v not in defrhs]
dead=[v for v in tree if tree[v] is None and v in defrhs]
print('live %d dead %d'%(len(live),len(dead)))

t0=time.time(); PIN={}
for L in live:
    n,side=leafnode[L]
    ws=[d[side] for d in OUT[n]]        # (coord1, coord2)
    _,r0=run({L:1})
    Cs=[]
    for w in ws:
        _,r1=run({L:1,w:1})
        cand=None
        for a in appearP[w]:
            i=ridx[a]; f0=r0[i]%p; sl=(r1[i]-r0[i])%p
            if f0 and sl:
                c=(-f0)*pow(sl,p-2,p)%p
                if cand is None: cand=c
                elif cand!=c: cand='CONFLICT'
        Cs.append(cand)
    PIN[L]=(ws,Cs)
print('pins extracted in %.0fs; leaves with a missing/conflicting pin: %s'%(
    time.time()-t0,[L for L,(w,c) in PIN.items() if any(x is None or x=='CONFLICT' for x in c)]))
pickle.dump({'PIN':PIN,'live':live,'dead':dead,'leafnode':leafnode},open('pins.pkl','wb'))
