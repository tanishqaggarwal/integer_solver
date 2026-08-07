"""PER-PLACEMENT SWEEP: put the defect at each cascade pin atom in turn and price the placement.
   For each: region R, inside-atoms S, balance deficit |R|-|S|, rank of the realizable
   zero-collateral knob image, and the EXACT integer optimum (max rows zeroable) -> failing."""
import ev, fast, json, time, itertools
from fast import St, csup
from chain import close_trace
from collections import defaultdict
from fractions import Fraction
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
atom_eqs=defaultdict(set)
for i,(m,sq,tl) in enumerate(ev.eq_terms):
    for c,a in tl: atom_eqs[a].add(i)
BITS=json.load(open('bits.json')); ALL=set(BITS['A']+BITS['B'])
CHAIN=[a for a,X in json.load(open('chain.json'))]
FREE=set(ev.F['free0'])
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
def solve_int(A,b,rows,n):
    M=[[A[j][i] for j in range(n)]+[-b[i]] for i in rows]
    piv=[]; r=0
    for c in range(n):
        k=None
        for i in range(r,len(M)):
            if M[i][c]: k=i;break
        if k is None: continue
        M[r],M[k]=M[k],M[r]
        for i in range(len(M)):
            if i!=r and M[i][c]:
                f1=M[r][c]; f2=M[i][c]
                M[i]=[f1*M[i][j]-f2*M[r][j] for j in range(n+1)]
        piv.append(c); r+=1
        if r==len(M): break
    for i in range(r,len(M)):
        if all(M[i][j]==0 for j in range(n)) and M[i][n]!=0: return None
    t=[0]*n
    for i,c in enumerate(piv):
        num=M[i][n]; den=M[i][c]
        if den==0 or num%den: return None
        t[c]=num//den
    for i in rows:
        if b[i]+sum(t[j]*A[j][i] for j in range(n))!=0: return None
    return t
def price(st,tag):
    NZ=set(st.nz())
    if not NZ: return None
    R=set()
    for a in NZ: R|=atom_eqs[a]
    S=sorted(a for a in set().union(*[set(q for c,q in ev.eq_terms[e][2]) for e in R]) if atom_eqs[a]<=R)
    Sset=set(S)
    cands=set()
    for a in S: cands.update([ev.F['free0'][b] for b in _bits(csup.get(a,0))] if a in csup else [])
    knobs=[]
    for X in sorted(cands):
        if X not in FREE: continue
        g=st.clone().set_free({X:st.fv.get(X,0)+1})
        d={a:g.av[a]-st.av[a] for a in g.av if g.av[a]!=st.av[a]}
        if d and all(a in Sset for a in d): knobs.append(X)
    Rl=sorted(R); b=[inner(st,e) for e in Rl]
    A=[]
    for X in knobs:
        g=st.clone().set_free({X:st.fv.get(X,0)+1})
        A.append([inner(g,e)-b[i] for i,e in enumerate(Rl)])
    n=len(knobs)
    rk=rank_q([[A[j][i] for j in range(n)] for i in range(len(Rl))],n) if n else 0
    z0=sum(1 for x in b if x==0); best=z0
    if n:
        for size in range(min(n,8),0,-1):
            if size<=best: break
            for rows in itertools.combinations(range(len(Rl)),size):
                t=solve_int(A,b,list(rows),n)
                if t is None: continue
                z=sum(1 for i in range(len(Rl)) if b[i]+sum(t[j]*A[j][i] for j in range(n))==0)
                if z>best: best=z
    outside=len([e for e in st.fails if e not in R])
    return dict(tag=tag,nz=len(NZ),R=len(Rl),S=len(S),deficit=len(Rl)-len(S),knobs=n,rank=rk,
                zeroable=best,failing=len(Rl)-best+outside,outside=outside,score=39033-(len(Rl)-best+outside))
def _bits(x):
    o=[]
    while x:
        q=x&-x; o.append(q.bit_length()-1); x^=q
    return o
if __name__=='__main__':
    st0=St({}); bsel=BITS['A'][0]
    rows=[]
    t0=time.time()
    for P in CHAIN:
        st=st0.clone().set_free({bsel:1})
        out,ok,tr,fr=close_trace(st,frozen=set(ALL),skip=frozenset([P]))
        r=price(out,'a%d'%P)
        if r: rows.append(r); print('  a%-6d nz=%-2d |R|=%-3d |S|=%-3d deficit=%-3d knobs=%-3d rank=%-3d zeroable=%-3d failing=%-3d score=%d  %.0fs'%(
            P,r['nz'],r['R'],r['S'],r['deficit'],r['knobs'],r['rank'],r['zeroable'],r['failing'],r['score'],time.time()-t0),flush=True)
    rows.sort(key=lambda d:-d['score'])
    print('\n=== PER-PLACEMENT TABLE (best first) ===')
    print('%-10s %4s %5s %5s %8s %6s %5s %8s %8s %7s'%('carrier','nz','|R|','|S|','deficit','knobs','rank','zeroable','failing','score'))
    for r in rows:
        print('%-10s %4d %5d %5d %8d %6d %5d %8d %8d %7d'%(r['tag'],r['nz'],r['R'],r['S'],r['deficit'],r['knobs'],r['rank'],r['zeroable'],r['failing'],r['score']))
    json.dump(rows,open('sweep.json','w'),indent=1)
    print('\nbest placement in the cascade: %d   deliverable: 39026'%rows[0]['score'])
