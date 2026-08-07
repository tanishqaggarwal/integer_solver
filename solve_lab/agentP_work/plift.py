#!/usr/bin/env python3
"""Agent P: integer lift constructor.
Walks the SLP in topological order; each atom is solved for its single unknown
variable over Z (recording any divisibility failure).  Reports exactly where the
integer lift breaks, if anywhere."""
import pickle,sys,json
from collections import defaultdict,Counter
sys.set_int_max_str_digits(10**7)
W='/home/user/integer_solver/solve_lab/agentP_work/'
D=pickle.load(open(W+'model4.pkl','rb')); AP=D['AP']
S=pickle.load(open(W+'slp.pkl','rb')); topo=S['topo']; outof=S['outof']
NV=38748
P=115792089237316195423570985008687907853269984665640564039457584007908834671663

def build(selset, verbose=True):
    """selset: set of selector variables to switch ON.  Returns (val, failures)."""
    val=[None]*NV
    import pfold as F
    for s in F.SEL: val[s]= 1 if s in selset else 0
    fails=[]; divfail=[]
    pend=list(range(len(topo)))
    for rnd in range(4):
        nxt=[]
        for i in pend:
            a=topo[i]; ap=AP[a]
            unk=sorted({x for m in ap for x in m if val[x] is None})
            if len(unk)==0:
                tot=0
                for m,c in ap.items():
                    t=c
                    for x in m: t*=val[x]
                    tot+=t
                if tot!=0: fails.append((i,a,tot))
                continue
            if len(unk)>1: nxt.append(i); continue
            y=unk[0]
            A=0; B=0
            for m,c in ap.items():
                k=list(m).count(y)
                if k==0:
                    t=c
                    for x in m: t*=val[x]
                    B+=t
                elif k==1:
                    t=c
                    for x in m:
                        if x!=y: t*=val[x]
                    A+=t
                else:                       # y^2 term -> boolean/idempotency atom
                    A=None; break
            if A is None:
                val[y]=0 if B==0 else None
                if val[y] is None: nxt.append(i)
                continue
            if A==0:
                if B!=0: fails.append((i,a,B))
                val[y]=0
                continue
            if (-B)%A!=0:
                divfail.append((i,a,A,B)); val[y]=(-B)//A       # nearest; atom will be nonzero
            else:
                val[y]=(-B)//A
        pend=nxt
        if not pend: break
    if verbose:
        print('  unresolved atoms after 4 passes:',len(pend))
        print('  divisibility failures (integer lift obstructions):',len(divfail))
        print('  atoms evaluating nonzero:',len(fails))
    for x in range(NV):
        if val[x] is None: val[x]=0
    return val,fails,divfail,pend

def score(val):
    D2=pickle.load(open(W+'model4.pkl','rb')); rows=D2['rows']
    av=[0]*len(AP)
    for i,ap in enumerate(AP):
        t=0
        for m,c in ap.items():
            u=c
            for x in m: u*=val[x]
            t+=u
        av[i]=t
    bad=[]
    for ei,r in enumerate(rows):
        L=sum(c*av[a] for c,a in r['row'])
        if L!=0: bad.append(ei)
    return bad,av

if __name__=='__main__':
    import pfold as F
    print("=== configuration: ALL SELECTORS OFF ===")
    val,fails,divfail,pend=build(set())
    bad,av=score(val)
    print('  equations failing:',len(bad),'-> score %d/39033'%(39033-len(bad)))
    print('  nonzero atoms:',sum(1 for x in av if x),
          ' at SLP pos:',sorted(i for i,a in enumerate(topo) if av[a])[:12])
    json.dump({'x%d'%i:str(v) for i,v in enumerate(val)},open(W+'lift_allzero.json','w'))
