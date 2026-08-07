"""HANDLE-CARRIER SWEEP.  A handle carries the defect in multiples of p: break its carrier atom by
   moving the handle one step (delta = m*p), freeze the handle, re-close everything else, and price.
   Stage A (cheap, ALL carriers): nz, |R|, |S|, deficit, knobs, rank.
   Stage B (expensive): exact integer optimum, only where rank is close to or above the deficit."""
import ev, fast, json, time, itertools, sys
from fast import St, csup
from chain import close_trace
from collections import defaultdict
from fractions import Fraction
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
atom_eqs=defaultdict(set)
for i,(m,sq,tl) in enumerate(ev.eq_terms):
    for c,a in tl: atom_eqs[a].add(i)
BITS=json.load(open('bits.json')); ALL=set(BITS['A']+BITS['B'])
FREE=set(ev.F['free0'])
H=[(int(X),int(a)) for X,a,m in json.load(open('handles.json'))]
NAMED=[(7497,688),(11436,1618),(22820,30982),(14393,30980)]   # non-solo, priced separately
def _bits(x):
    o=[]
    while x:
        q=x&-x; o.append(q.bit_length()-1); x^=q
    return o
FR0=ev.F['free0']
def inner(st,e):
    m,sq,tl=ev.eq_terms[e]; t=0
    for c,a in tl:
        x=st.av.get(a)
        if x: t+=c*x
    return t
def rank_q(rows,n):
    rr=[[Fraction(x) for x in r] for r in rows]; piv=0
    for c in range(n):
        k=None
        for i in range(piv,len(rr)):
            if rr[i][c]!=0: k=i;break
        if k is None: continue
        rr[piv],rr[k]=rr[k],rr[piv]; pv=rr[piv][c]
        for i in range(len(rr)):
            if i!=piv and rr[i][c]!=0:
                f=rr[i][c]/pv; rr[i]=[rr[i][j]-f*rr[piv][j] for j in range(n)]
        piv+=1
        if piv==len(rr): break
    return piv
# base: the closed one-selector state (39,018)
st0=St({}); bsel=BITS['A'][0]
base,ok,tr,frz=close_trace(st0.clone().set_free({bsel:1}),frozen=set(ALL))
print('base closed state score',base.score(),flush=True)
def stageA(X,a):
    g=base.clone().set_free({X:base.fv.get(X,0)+1})
    out,ok2,tr2,_=close_trace(g,frozen=set(ALL)|{X})
    NZ=set(out.nz())
    if not NZ: return None
    R=set()
    for q in NZ: R|=atom_eqs[q]
    Rl=sorted(R)
    S=sorted(q for q in set().union(*[set(z for c,z in ev.eq_terms[e][2]) for e in Rl]) if atom_eqs[q]<=R)
    Sset=set(S)
    cands=set()
    for q in S:
        if q in csup: cands.update(FR0[b] for b in _bits(csup[q]))
    knobs=[]
    for Y in sorted(cands):
        if Y not in FREE: continue
        h=out.clone().set_free({Y:out.fv.get(Y,0)+1})
        d={q:h.av[q]-out.av[q] for q in h.av if h.av[q]!=out.av[q]}
        if d and all(q in Sset for q in d): knobs.append(Y)
    b=[inner(out,e) for e in Rl]
    A=[]
    for Y in knobs:
        h=out.clone().set_free({Y:out.fv.get(Y,0)+1})
        A.append([inner(h,e)-b[i] for i,e in enumerate(Rl)])
    n=len(knobs)
    rk=rank_q([[A[j][i] for j in range(n)] for i in range(len(Rl))],n) if n else 0
    z0=sum(1 for x in b if x==0)
    outside=len([e for e in out.fails if e not in R])
    return dict(handle=X,carrier=a,nz=len(NZ),R=len(Rl),S=len(S),deficit=len(Rl)-len(S),
                knobs=n,rank=rk,z0=z0,failing=len(Rl)-z0+outside,score=39033-(len(Rl)-z0+outside),
                A=A,b=b,Rl=Rl,state=out)
res=[]; t0=time.time()
TARGETS=NAMED+H
for i,(X,a) in enumerate(TARGETS):
    try:
        r=stageA(X,a)
    except Exception as e:
        print('  x_%d ERROR %s'%(X,e),flush=True); continue
    if r is None: continue
    res.append({k:v for k,v in r.items() if k not in ('A','b','Rl','state')})
    if r['rank']>r['deficit'] or r['score']>39018:
        print('  *** x_%-6d carrier a%-6d rank=%d > deficit=%d  score=%d'%(X,a,r['rank'],r['deficit'],r['score']),flush=True)
    if i%50==0:
        print('  %d/%d  %.0fs  best score so far %d'%(i,len(TARGETS),time.time()-t0,max(z['score'] for z in res)),flush=True)
        json.dump(res,open('hsweep_partial.json','w'))
json.dump(res,open('hsweep.json','w'))
print('\nstage A complete: %d carriers priced in %.0fs'%(len(res),time.time()-t0))
win=[r for r in res if r['rank']>r['deficit']]
print('carriers with rank > deficit: %d'%len(win))
for r in sorted(res,key=lambda z:-z['score'])[:15]:
    print('  x_%-6d a%-6d nz=%-2d |R|=%-3d |S|=%-3d deficit=%-3d knobs=%-3d rank=%-3d failing=%-3d score=%d'%(
        r['handle'],r['carrier'],r['nz'],r['R'],r['S'],r['deficit'],r['knobs'],r['rank'],r['failing'],r['score']))
