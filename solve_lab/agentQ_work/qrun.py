#!/usr/bin/env python3
"""Q-10: run the REAL gate DAG for a chosen ON-set S and compare the root wire with fold(S).

Everything is done mod p (the circuit's integer lifting is irrelevant to the point values).
Known at the start: the 256 selector bits, and the 256 leaf coordinate wires (their pinned
constants).  Every other value is produced either by an actual gate of EQUATIONS.txt or by a
chord-law gadget, whose two output wires (u3,y3) are FREE inputs determined by its own residuals:
      lam = dy/dx,   u3 = lam^2 - ua - ub - K,   y3 = lam*(ub - u3) - yb
No group-theory input is used in the propagation -- the group law is only used afterwards, to
compute fold(S) independently, and the two are compared.
"""
import json,collections,random,sys,os,re
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from qgrp import add,mul,oncur,p,cs
K=97553848499418123410591666447050222001188385549510401465815187079080512838891
ROOTX,ROOTY=24468,18956
G=[json.loads(l) for l in open('atoms/gates.jsonl')]
alld=collections.defaultdict(list)
for d in G: alld[d['t']].append(d)
VAR=re.compile(r'x_(\d+)')
comp={}
for t,ds in alld.items():
    comp[t]=[(compile(VAR.sub(r'V[\1]',d['rhs']),'<g>','eval'),tuple(sorted(set(d['vids'])))) for d in ds]
uses=collections.defaultdict(set)
for t,ds in alld.items():
    for d in ds:
        for v in d['vids']: uses[v].add(t)
leaf={int(g):v for g,v in json.load(open('qleaf.json')).items()}
lad=json.load(open('qladder.json'))
e2s={int(k):v for k,v in lad['exp2sel'].items()}
LP={int(g):(int(v[0]),int(v[1])) for g,v in leaf.items()}      # curve point
LW={int(g):(v[2],v[3]) for g,v in leaf.items()}                # (Xwire, Ywire)
ST=[s for s in json.load(open('qstages.json'))['stages'] if 'u3' in s]
gad_by_in=collections.defaultdict(list)
for i,s in enumerate(ST):
    for w in (s['ua'],s['ub'],s['ya'],s['yb']): gad_by_in[w].append(i)
NV=38748
def run(S, leafmode='pin'):
    V=[None]*NV
    for g in LP:
        V[g]=1 if g in S else 0
        wx,wy=LW[g]
        if leafmode=='pin' or g in S:
            V[wx]=(LP[g][0]-cs)%p; V[wy]=LP[g][1]%p
        else:
            V[wx]=0; V[wy]=0
    ns0={'V':V,'__builtins__':{}}
    for t,cs_ in comp.items():                 # constant gates (rhs has no variable)
        for c,vids in cs_:
            if not vids and V[t] is None: V[t]=eval(c,ns0)%p; break
    q=collections.deque(w for w in range(NV) if V[w] is not None)
    done=set(); gdone=set()
    ns={'V':V,'__builtins__':{}}
    while q:
        w=q.popleft()
        for t in uses[w]:
            if V[t] is not None: continue
            for c,vids in comp[t]:
                if all(V[v] is not None for v in vids):
                    V[t]=eval(c,ns)%p; q.append(t); break
        for gi in gad_by_in[w]:
            if gi in gdone: continue
            s=ST[gi]
            if any(V[s[k]] is None for k in ('ua','ub','ya','yb')): continue
            if V[s['u3']] is not None: gdone.add(gi); continue
            dx=(V[s['ua']]-V[s['ub']])%p; dy=(V[s['ya']]-V[s['yb']])%p
            if dx==0: gdone.add(gi); continue
            lam=dy*pow(dx,p-2,p)%p
            u3=(lam*lam-V[s['ua']]-V[s['ub']]-K)%p
            y3=(lam*(V[s['ub']]-u3)-V[s['yb']])%p
            V[s['u3']]=u3; V[s['y3']]=y3; gdone.add(gi)
            q.append(s['u3']); q.append(s['y3'])
    return V,len(gdone)
def fold(S):
    F=None
    for g in S: F=add(F,LP[g])
    return F
random.seed(3)
sel=sorted(LP)
print('%-6s %-8s %-10s %-8s %s'%('|S|','gadgets','rootX set','matches','note'))
for n in (1,2,5,17,64,128,200,256):
    S=set(random.sample(sel,n))
    V,ng=run(S)
    F=fold(S)
    got=V[ROOTX]
    want=(F[0]-cs)%p if F else None
    ok = got is not None and want is not None and got==want and V[ROOTY]==F[1]%p
    print('%-6d %-8d %-10s %-8s %s'%(n,ng,got is not None,ok,'' if ok else 'MISMATCH/UNSET'))
