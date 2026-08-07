"""TWO-SELECTOR CLOSURE in my own frame, and stage-A pricing of the resulting region.
   Question: does the zero-collateral knob-image RANK rise faster than the balance DEFICIT when a
   second same-tree selector is on?  Criterion unchanged: rank > deficit."""
import ev, fast, json, time, itertools, sys
from fast import St, csup
from chain import close_trace
from collections import defaultdict, Counter
from fractions import Fraction
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
atom_eqs=defaultdict(set)
for i,(m,sq,tl) in enumerate(ev.eq_terms):
    for c,a in tl: atom_eqs[a].add(i)
BITS=json.load(open('bits.json')); UW=json.load(open('uw.json'))
ALL=set(BITS['A']+BITS['B']); FREE=set(ev.F['free0']); FR0=ev.F['free0']
def _bits(x):
    o=[]
    while x:
        q=x&-x; o.append(q.bit_length()-1); x^=q
    return o
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
def price(out,tag):
    NZ=set(out.nz())
    if not NZ: return dict(tag=tag,nz=0,R=0,S=0,deficit=0,knobs=0,rank=0,failing=len(out.fails),score=out.score())
    R=set()
    for q in NZ: R|=atom_eqs[q]
    Rl=sorted(R)
    S=sorted(q for q in set().union(*[set(z for c,z in ev.eq_terms[e][2]) for e in Rl]) if atom_eqs[q]<=R)
    Sset=set(S); cands=set()
    for q in S:
        if q in csup: cands.update(FR0[b] for b in _bits(csup[q]))
    knobs=[]
    for Y in sorted(cands):
        if Y not in FREE: continue
        h=out.clone().set_free({Y:out.fv.get(Y,0)+1})
        d={q:h.av[q]-out.av[q] for q in h.av if h.av[q]!=out.av[q]}
        if d and all(q in Sset for q in d): knobs.append(Y)
    b=[inner(out,e) for e in Rl]; A=[]
    for Y in knobs:
        h=out.clone().set_free({Y:out.fv.get(Y,0)+1})
        A.append([inner(h,e)-b[i] for i,e in enumerate(Rl)])
    n=len(knobs)
    rk=rank_q([[A[j][i] for j in range(n)] for i in range(len(Rl))],n) if n else 0
    z0=sum(1 for x in b if x==0)
    outside=len([e for e in out.fails if e not in R])
    return dict(tag=tag,nz=len(NZ),R=len(Rl),S=len(S),deficit=len(Rl)-len(S),knobs=n,rank=rk,
                z0=z0,failing=len(Rl)-z0+outside,score=39033-(len(Rl)-z0+outside),knoblist=knobs)
st0=St({}); U=UW['U']; W=UW['W']
# 1-selector reference
r1=price(close_trace(st0.clone().set_free({U[0]:1}),frozen=set(ALL))[0],'ONE u[0]')
print('ONE-SELECTOR  : nz=%d |R|=%d |S|=%d deficit=%d knobs=%d rank=%d failing=%d score=%d'%(
    r1['nz'],r1['R'],r1['S'],r1['deficit'],r1['knobs'],r1['rank'],r1['failing'],r1['score']),flush=True)
print('              knobs:',r1['knoblist'],flush=True)
rows=[r1]
t0=time.time()
# 2 same-tree selectors, varied LCA depth: adjacent in the u list, far apart, and cross u/w
pairs=[('same-u adjacent',U[0],U[1]),('same-u near',U[0],U[3]),('same-u mid',U[0],U[20]),
       ('same-u far',U[0],U[60]),('same-u other',U[10],U[50]),('same-u other2',U[30],U[70]),
       ('same-w adjacent',W[0],W[1]),('same-w far',W[0],W[60]),
       ('cross u/w',U[0],W[0]),('cross u/w 2',U[13],W[21])]
for tag,b1,b2 in pairs:
    out,ok,tr,_=close_trace(st0.clone().set_free({b1:1,b2:1}),frozen=set(ALL))
    r=price(out,'%s (%d,%d)'%(tag,b1,b2)); rows.append(r)
    print('%-22s: nz=%-2d |R|=%-3d |S|=%-3d deficit=%-3d knobs=%-3d rank=%-3d failing=%-3d score=%d  gap=%d  %.0fs'%(
        tag,r['nz'],r['R'],r['S'],r['deficit'],r['knobs'],r['rank'],r['failing'],r['score'],
        r['deficit']-r['rank'],time.time()-t0),flush=True)
    if r['rank']>r['deficit']: print('   *** RANK > DEFICIT ***',flush=True)
    if r['score']>39018: print('   *** beats the one-selector 39,018 ***',flush=True)
json.dump([{k:v for k,v in r.items() if k!='knoblist'} for r in rows],open('two.json','w'),indent=1)
best=max(rows,key=lambda z:z['score'])
print('\nbest two-selector score: %d (one-selector: %d, deliverable: 39026)'%(best['score'],r1['score']))
