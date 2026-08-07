#!/usr/bin/env python3
"""Agent P: the audit checks."""
import pickle,sys,json,random
from collections import Counter,defaultdict
sys.set_int_max_str_digits(10**7)
sys.setrecursionlimit(20000)
import pfold as F
P,Q=F.P,F.Q
W=F.W
g=[0]*38748
for k,v in json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')).items(): g[int(k[2:])]=int(v)

print("="*70)
print("A. STAGE COUNT / SHAPE")
kinds=Counter()
for r in F.SRC:
    kinds[tuple(sorted(x[0] for x in r))]+=1
print("  law-blocks total (Q-gates):",383,"  with mux:",len(F.stages),"  + 1 root")
print("  stage input-kind profile:",dict(kinds))
real=[j for j,r in enumerate(F.SRC) if all(x[0]!='0' for x in r)]
print("  stages with BOTH inputs non-identity:",len(real))
dead=[j for j,s in enumerate(F.stages) if len(F.supp[j])==0]
print("  dead stages (empty leaf support):",len(dead))
print("  live stages:",382-len(dead))

# contracted tree: skip identity stages
def contract():
    ch={}
    def node(j):
        kids=[]
        for k in F.SRC[j]:
            if k[0]=='S':
                if len(F.supp[k[1]])==0: continue
                kids.append(('S',k[1]))
            elif k[0]=='L': kids.append(('L',k[1]))
        return kids
    return {j:node(j) for j in range(382)}
CH=contract()
def cdepth(j,memo={}):
    if j in memo: return memo[j]
    kids=CH[j]
    if not kids: d=0
    else:
        ds=[]
        for k in kids:
            if k[0]=='S':
                ds.append(cdepth(k[1]) + (1 if len(CH[k[1]])>1 else 0))
            else: ds.append(0)
        d=max(ds)+(1 if len(kids)>1 else 0)
    memo[j]=d; return d
merges=[j for j in range(382) if len(CH[j])>1]
print("  MERGE stages (2 live inputs):",len(merges))
print("  pass-through stages (1 live input):",sum(1 for j in range(382) if len(CH[j])==1))
rootdepth=max(cdepth(380),cdepth(381))+1
print("  contracted merge-depth (root incl.):",rootdepth)
print("  raw depth (all stages):",max(F.depth_ if False else 0,0) if False else max(pickle.load(open(W+'topo.pkl','rb'))['depth'])+1)

print("="*70)
print("B. ROOT LAW verified against the deliverable's own values")
X=(g[F.stages[381]['out'][0]]%P, g[F.stages[381]['out'][1]]%P)
Y=(g[F.stages[380]['out'][0]]%P, g[F.stages[380]['out'][1]]%P)
r=F.law(X,Y); r2=F.law(Y,X)
tgt=(g[13682]%P,g[37892]%P)
print("  law(S381,S380) == deliverable root pair:",r==tgt)
print("  law(S380,S381) == deliverable root pair:",r2==tgt)
print("  deliverable root pair == TARGET:",tgt==(F.TGT_X,F.TGT_Y))

print("="*70)
print("C. LAW INVERTIBILITY on random triples")
rnd=random.Random(20260807)
okf=okb=0
for _ in range(200):
    x1=rnd.randrange(P); y1=rnd.randrange(P); x2=rnd.randrange(P); y2=rnd.randrange(P)
    Z=F.law((x1,y1),(x2,y2))
    if Z is None: continue
    # invert: given X=(x1,y1) and Z, recover Y
    x3,y3=Z
    # lam = (y3+y1)/(x1-x3)  from y3 = lam(x1-x3)-y1
    d=(x1-x3)%P
    if d==0: continue
    lam=((y3+y1)%P)*pow(d,P-2,P)%P
    x2r=(lam*lam-x1-x3-Q)%P
    y2r=(lam*(x1-x2r)-y1)%P
    if (x2r,y2r)==(x2,y2): okb+=1
    if F.law((x1,y1),(x2r,y2r))==Z: okf+=1
print("  inverse recovers Y exactly:",okb,"/200 ; forward re-check:",okf)

print("="*70)
print("D. HOW MANY BOOLEAN FREEDOMS?")
D=pickle.load(open(W+'model4.pkl','rb')); AP=D['AP']
bl=set()
for ap in AP:
    ks=list(ap.items())
    if len(ks)==2:
        (m1,c1),(m2,c2)=sorted(ks,key=lambda z:len(z[0]))
        if len(m1)==1 and len(m2)==2 and m2[0]==m2[1]==m1[0] and c1==-c2: bl.add(m1[0])
print("  boolean-constrained variables:",len(bl))
print("  of which are leaf selectors:",len(bl & set(F.SEL)))
print("  other booleans:",len(bl-set(F.SEL)))
