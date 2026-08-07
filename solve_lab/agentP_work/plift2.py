#!/usr/bin/env python3
"""Agent P: integer lift constructor, worklist fixed point.
Seeds = 256 selectors + 512 leaf coords + 764 block law outputs.  Everything else is
propagated by solving each atom for its single remaining unknown OVER Z.  Any place the
solve needs a division that does not come out exactly is recorded as a lift obstruction."""
import pickle,sys,json
from collections import defaultdict,Counter,deque
sys.set_int_max_str_digits(10**7)
W='/home/user/integer_solver/solve_lab/agentP_work/'
P=115792089237316195423570985008687907853269984665640564039457584007908834671663
D=pickle.load(open(W+'model4.pkl','rb')); AP=D['AP']; rows=D['rows']
S=pickle.load(open(W+'slp.pkl','rb')); topo=S['topo']; outof=S['outof']
B=pickle.load(open(W+'blocks.pkl','rb'))
LEAVES=pickle.load(open(W+'leaves.pkl','rb'))
NV=38748
import pfold as F
v2a=defaultdict(list)
for i,ap in enumerate(AP):
    for x in set(y for m in ap for y in m): v2a[x].append(i)
LK=defaultdict(list)
for a,b,k in LEAVES: LK[a].append((b,k))

def build(selset):
    val=[None]*NV
    for s in F.SEL: val[s]=1 if s in selset else 0
    for s in F.SEL:
        on = 1 if s in selset else 0
        for (b,k) in LK[s]: val[b]= k if on else 0
    # block law outputs: honest fold value, lifted into [0,P)
    memo={}
    def ev(j):
        if j in memo: return memo[j]
        vals=[]
        for k in F.SRC[j]:
            if k[0]=='S': vals.append(ev(k[1]))
            elif k[0]=='L':
                bit=1 if k[1] in selset else 0
                vals.append((k[2],1) if bit else ((0,0),0))
            else: vals.append(((0,0),0))
        (X,a),(Y,b)=vals
        if a and b:
            r=F.law(X,Y); res=(r,1) if r else (None,1)
        elif a: res=(X,1)
        elif b: res=(Y,1)
        else: res=((0,0),0)
        memo[j]=res; return res
    sys.setrecursionlimit(20000)
    for j in range(len(B)): ev(j)
    for j,b in enumerate(B):
        z,live=memo[j]
        live2 = 1 if all(k[0]!='0' and (k[0]!='L' or k[1] in selset) for k in F.SRC[j]) else 0
        if z is not None and live2:
            val[b['i5']],val[b['i6']]=z[0],z[1]
        else:
            val[b['i5']]=0; val[b['i6']]=0
    obstruct=[]; nonzero=[]
    q=deque(range(len(AP)))
    inq=[True]*len(AP)
    while q:
        i=q.popleft(); inq[i]=False; ap=AP[i]
        unk=sorted({x for m in ap for x in m if val[x] is None})
        if len(unk)!=1: continue
        y=unk[0]
        A=0;B_=0;sq=False
        for m,c in ap.items():
            k=list(m).count(y)
            t=c
            if k==0:
                for x in m: t*=val[x]
                B_+=t
            elif k==1:
                for x in m:
                    if x!=y: t*=val[x]
                A+=t
            else: sq=True
        if sq:
            if A==0 and B_==0: val[y]=0
            else: continue
        elif A==0:
            continue
        else:
            if (-B_)%A!=0:
                obstruct.append((i,A,B_)); val[y]=(-B_)//A
            else: val[y]=(-B_)//A
        for j2 in v2a[y]:
            if not inq[j2]: inq[j2]=True; q.append(j2)
    und=sum(1 for x in range(NV) if val[x] is None)
    for x in range(NV):
        if val[x] is None: val[x]=0
    return val,obstruct,und

def score(val):
    av=[0]*len(AP)
    for i,ap in enumerate(AP):
        t=0
        for m,c in ap.items():
            u=c
            for x in m: u*=val[x]
            t+=u
        av[i]=t
    bad=[ei for ei,r in enumerate(rows) if sum(c*av[a] for c,a in r['row'])!=0]
    return bad,av

if __name__=='__main__':
    for name,sel in [('ALL OFF',set()),('one leaf ON',{F.SEL[0]})]:
        print('=== configuration: %s ==='%name)
        val,obs,und=build(sel)
        bad,av=score(val)
        pos={a:i for i,a in enumerate(topo)}
        nz=sorted(pos[a] for a in range(len(AP)) if av[a])
        print('  undetermined vars:',und,'  integer-division obstructions:',len(obs))
        print('  nonzero atoms:',len(nz),nz[:14])
        print('  equations failing:',len(bad),'  SCORE %d/39033'%(39033-len(bad)))
        json.dump({'x%d'%i:str(v) for i,v in enumerate(val)},
                  open(W+'lift_%s.json'%name.replace(' ','_'),'w'))
