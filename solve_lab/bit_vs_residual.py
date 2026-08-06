import json
import heal_harness as H
from collections import defaultdict
p=H.p
pinrec=json.load(open('pinrec.json'))
bysel=defaultdict(list)
for i,sel,tgt,const,coef,handle in pinrec: bysel[sel].append((tgt,const%p,handle))
bits=sorted(bysel)
d=H.loadd('g1g2_closed.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
H.forward()
V=H.val
atoms=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f: atoms.append(json.loads(line))
RES=[7450,7452,44342,45677]
def rvec():
    out=[]
    for ai in RES:
        s=0
        for vs,c in atoms[ai]['poly']:
            m=c
            for vi in vs: m*=V[vi]
            s+=m
        out.append(s%p)
    return out
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':V,'__builtins__':{}}
def fwd_from(knobs):
    aff=set()
    for w in knobs: aff|=set(desc_of[w])
    for k in sorted(aff): V[H.order[k]]=eval(H.gcode[k],ns)
base=rvec()
print(f"current residual vector (mod p): {[str(x)[:12]+'..' for x in base]}")
# measure each bit's effect on the 4 residuals (toggle 0->1, activate loads)
effects={}
for b in bits:
    saved={b:V[b]}; V[b]=1-V[b]
    for tgt,const,handle in bysel[b]:
        saved[tgt]=V[tgt]; V[tgt]=const if V[b]==1 else 0
    knobs=[b]+[t for t,_,_ in bysel[b]]
    fwd_from(knobs)
    e=[(rvec()[j]-base[j])%p for j in range(4)]
    for v,x in saved.items(): V[v]=x
    fwd_from(knobs)
    if any(e): effects[b]=e
print(f"bits affecting the 4 residuals: {len(effects)}")
# build matrix over GF(p), check if -base is in column span
import itertools
bl=list(effects)
if bl:
    # Gaussian elimination mod p to see rank and whether target reachable
    # rows=4 residuals, cols=bits. Solve E * x = -base mod p
    rows=4; cols=len(bl)
    A=[[effects[bl[j]][i] for j in range(cols)] for i in range(rows)]
    tgt=[(-base[i])%p for i in range(rows)]
    # augmented gaussian
    M=[A[i][:]+[tgt[i]] for i in range(rows)]
    pr=0
    for col in range(cols):
        piv=None
        for r in range(pr,rows):
            if M[r][col]%p!=0: piv=r;break
        if piv is None: continue
        M[pr],M[piv]=M[piv],M[pr]
        inv=pow(M[pr][col],-1,p)
        M[pr]=[(x*inv)%p for x in M[pr]]
        for r in range(rows):
            if r!=pr and M[r][col]%p!=0:
                f=M[r][col]; M[r]=[(M[r][k]-f*M[pr][k])%p for k in range(cols+1)]
        pr+=1
    # check consistency: any row 0...0 | nonzero
    incons=False
    for r in range(rows):
        if all(M[r][k]%p==0 for k in range(cols)) and M[r][cols]%p!=0: incons=True
    print(f"rank={pr}, target reachable (consistent): {not incons}")
    if not incons:
        print("*** RESIDUAL IS IN BIT-SPAN — a bit-combination zeros it mod p! ***")
