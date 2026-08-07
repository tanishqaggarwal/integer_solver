#!/usr/bin/env python3
"""Agent P: integer lift constructor, seeding fix.

Fix over plift2: when an atom holds exactly TWO unknowns and one of them is a handle
variable h (a variable defined by an atom `h - (P-alias)*u`), set h = 0.  That is the
general form of "seed every mod-P copy target to equal its source exactly over Z":
h = 0 forces u = 0 in h's own definition atom and makes the copy exact, unblocking
the cascade.  Rule 1 (single unknown) always fires first, so handles whose residual R
is genuinely determined are still solved as h = R/c and their divisibility is recorded.
"""
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

# --- handle set (same definition as check-in 11) ---
par=list(range(NV))
def _f(a):
    while par[a]!=a: par[a]=par[par[a]]; a=par[a]
    return a
pinP=None
for ap in AP:
    it=list(ap.items())
    if len(it)==2:
        (m1,c1),(m2,c2)=sorted(it,key=lambda z:len(z[0]))
        if len(m1)==1 and len(m2)==1 and c1==-c2:
            a,b=_f(m1[0]),_f(m2[0])
            if a!=b: par[a]=b
        elif m1==() and len(m2)==1 and abs(c1)==P: pinP=m2[0]
PV={x for x in range(NV) if _f(x)==_f(pinP)}
HDEF={}
for i,ap in enumerate(AP):
    o=outof[i]
    if o is None: continue
    q=[m for m in ap if len(m)==2]
    if len(q)==1 and len(ap)==2 and (q[0][0] in PV)!=(q[0][1] in PV):
        HDEF[o]=i
HANDLE=set(HDEF)

def build(selset):
    val=[None]*NV
    for s in F.SEL: val[s]=1 if s in selset else 0
    for s in F.SEL:
        on=1 if s in selset else 0
        for (b,k) in LK[s]: val[b]= k if on else 0
    memo={}
    def ev(j):
        if j in memo: return memo[j]
        vals=[]
        for k in F.SRC[j]:
            if k[0]=='S': vals.append(ev(k[1]))
            elif k[0]=='L':
                vals.append((k[2],1) if k[1] in selset else ((0,0),0))
            else: vals.append(((0,0),0))
        (X,a),(Y,b)=vals
        if a and b:
            r=F.law(X,Y); res=(r,1) if r else (None,1)
        elif a: res=(X,1)
        elif b: res=(Y,1)
        else: res=((0,0),0)
        memo[j]=res; return res
    sys.setrecursionlimit(20000)
    live={}
    for j in range(len(B)):
        z,l=ev(j)
        lv = all((k[0]=='S' and len(F.supp[k[1]] if False else [])>=0 and ev(k[1])[1]==1) or
                 (k[0]=='L' and k[1] in selset) for k in F.SRC[j])
        live[j]=1 if lv else 0
        b=B[j]
        if lv and z is not None: val[b['i5']],val[b['i6']]=z[0],z[1]
        else: val[b['i5']],val[b['i6']]=0,0
    obstruct=[]; setzero=0
    q=deque(range(len(AP))); inq=[True]*len(AP)
    while q:
        i=q.popleft(); inq[i]=False; ap=AP[i]
        unk=sorted({x for m in ap for x in m if val[x] is None})
        if len(unk)==2:
            hs=[y for y in unk if y in HANDLE]
            if len(hs)==1:
                val[hs[0]]=0; setzero+=1
                for j2 in v2a[hs[0]]:
                    if not inq[j2]: inq[j2]=True; q.append(j2)
                if not inq[i]: inq[i]=True; q.append(i)
            continue
        if len(unk)!=1: continue
        y=unk[0]; A=0; Bc=0; sq=False
        for m,c in ap.items():
            k=list(m).count(y); t=c
            if k==0:
                for x in m: t*=val[x]
                Bc+=t
            elif k==1:
                for x in m:
                    if x!=y: t*=val[x]
                A+=t
            else: sq=True
        if sq:
            if A==0 and Bc==0: val[y]=0
            else: continue
        elif A==0: continue
        else:
            if (-Bc)%A!=0: obstruct.append((i,y,A,Bc))
            val[y]=(-Bc)//A
        for j2 in v2a[y]:
            if not inq[j2]: inq[j2]=True; q.append(j2)
    und=[x for x in range(NV) if val[x] is None]
    for x in und: val[x]=0
    return val,obstruct,und,setzero,live

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
    pos={a:i for i,a in enumerate(topo)}
    # a configuration with a genuinely live merge: block 2 has two leaf inputs
    twoleaf=set()
    for j,row in enumerate(F.SRC):
        if all(k[0]=='L' for k in row):
            twoleaf={row[0][1],row[1][1]}; jj=j; break
    cfgs=[('all_off',set()),('one_leaf',{F.SEL[0]}),('two_leaf_block%d'%jj,twoleaf)]
    for name,sel in cfgs:
        val,obs,und,sz,live=build(sel)
        bad,av=score(val)
        nz=sorted(pos[a] for a in range(len(AP)) if av[a])
        print('=== %s ===' % name)
        print('  live blocks:',sum(live.values()),' handles forced to 0:',sz)
        print('  undetermined vars:',len(und),'  integer-division obstructions:',len(obs))
        print('  nonzero atoms:',len(nz),nz[:10])
        print('  equations failing:',len(bad),'  (score %d/39033)'%(39033-len(bad)))
        json.dump({'x%d'%i:str(v) for i,v in enumerate(val)},open(W+'lift3_%s.json'%name,'w'))
