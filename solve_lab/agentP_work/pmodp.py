#!/usr/bin/env python3
"""Agent P: mod-P propagator. Choose 256 selector bits -> solve the whole F_P system."""
import pickle,sys,json,random
from collections import defaultdict,Counter,deque
sys.set_int_max_str_digits(10**7)
W='/home/user/integer_solver/solve_lab/agentP_work/'
P=115792089237316195423570985008687907853269984665640564039457584007908834671663
NV=38748

_D=pickle.load(open(W+'model4.pkl','rb'))
AP=_D['AP']; ROWS=_D['rows']
_g=[0]*NV
for k,v in json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')).items(): _g[int(k[2:])]=int(v)
PVARS={x for x in range(NV) if _g[x]==P}
LEAVES=pickle.load(open(W+'leaves.pkl','rb'))          # (sel, coordvar, K)
SEL=sorted(set(t[0] for t in LEAVES))
LEAFOF=defaultdict(list)
for a,b,k in LEAVES: LEAFOF[a].append((b,k%P))

# reduce atoms mod P, drop P-vars
RED=[]
for ap in AP:
    d=defaultdict(int)
    for m,c in ap.items():
        if any(x in PVARS for x in m): continue
        d[m]=(d[m]+c)%P
    RED.append({m:c for m,c in d.items() if c})
var2at=defaultdict(list)
for i,r in enumerate(RED):
    s=set()
    for m in r: s.update(m)
    for x in s: var2at[x].append(i)

def solve(bits, verbose=True):
    """bits: dict selvar->0/1 (or list of 256). Returns val dict."""
    val={x:0 for x in PVARS}
    if isinstance(bits,(list,tuple)):
        bits={SEL[i]:b for i,b in enumerate(bits)}
    for s,b in bits.items(): val[s]=b%P
    q=deque(range(len(RED)))
    inq=[True]*len(RED)
    rounds=0
    while q:
        i=q.popleft(); inq[i]=False; rounds+=1
        r=RED[i]
        unk=set()
        for m in r:
            for x in m:
                if x not in val: unk.add(x)
        if len(unk)!=1: continue
        x=next(iter(unk))
        A=0;B=0; lin=True
        for m,c in r.items():
            k=m.count(x)
            if k==0:
                t=c
                for y in m: t=t*val[y]%P
                B=(B+t)%P
            elif k==1:
                t=c
                for y in m:
                    if y!=x: t=t*val[y]%P
                A=(A+t)%P
            else: lin=False; break
        if not lin: continue
        if A==0: continue
        val[x]=(-B)*pow(A,P-2,P)%P
        for j in var2at[x]:
            if not inq[j]: inq[j]=True; q.append(j)
    if verbose: print(f"  solved {len(val)}/{NV} vars, {rounds} pops")
    return val

def check(val):
    bad=[]
    for i,r in enumerate(RED):
        s=0; ok=True
        for m,c in r.items():
            t=c
            for x in m:
                if x not in val: ok=False;break
                t=t*val[x]%P
            if not ok: break
            s=(s+t)%P
        if ok and s: bad.append(i)
    return bad

if __name__=='__main__':
    print("selectors:",len(SEL))
    b0={s:(1 if _g[s] else 0) for s in SEL}
    print("deliverable ON:",sum(b0.values()))
    v=solve(b0)
    unset=[x for x in range(NV) if x not in v]
    print("unset vars:",len(unset))
    bad=check(v)
    print("violated congruences:",len(bad))
    mism=[x for x in v if _g[x]%P!=v[x]]
    print("mismatch vs deliverable mod P:",len(mism))
    pickle.dump({'unset':unset,'bad':bad},open(W+'modp_diag.pkl','wb'))
