"""Detach-set sweep.  ONE maximal frame with the whole pool detached; a subset D is emulated by
   re-attaching every pool variable NOT in D to its gate value.  Criterion: deficit < 4 or rank > 7."""
import ev, json, time, itertools, sys, ast, re
import frameB as FB
from frameB import Frame, State
from fractions import Fraction
from collections import defaultdict
VAR_RE=re.compile(r'x_(\d+)')
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
atom_eqs=defaultdict(set)
for i,(m,sq,tl) in enumerate(ev.eq_terms):
    for c,a in tl: atom_eqs[a].add(i)
POOL=json.load(open('pool.json'))
WIT=[642,28730,29854,31864]
ORDER=[v for v in ev.F['order'] if v in set(POOL)]          # topological order within the pool
DEFRHS={}
for v in POOL:
    a=ev.F['definer'][v]
    t=ast.parse(ev.atom_src[a],mode='eval').body
    DEFRHS[v]=compile(VAR_RE.sub(r'v[\1]',ast.unparse(t.right)),'<d>','eval')
t0=time.time()
fr=Frame(POOL)
print('maximal frame: free=%d checks=%d  built %.1fs'%(len(fr.free),len(fr.checks),time.time()-t0),flush=True)
W=json.load(open('../best/new_instance_partial_39026.json'))
wv=[0]*38748
for k,val in W.items(): wv[int(k[2:]) if k.startswith('x_') else int(k)]=int(val)
BASEFV={u:wv[u] for u in fr.free if wv[u]!=0}
def make(D):
    st=State(fr,BASEFV)
    Ds=set(D)
    for _ in range(3):
        ch={}
        for v in ORDER:
            if v in Ds: continue
            ns={'v':st.v,'__builtins__':{}}
            ch[v]=eval(DEFRHS[v],ns)
            st.set_free({v:ch[v]})
    return st
FREE=set(fr.free); FR0=fr.free
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
def price(st,tag):
    NZ=set(st.nz())
    if not NZ: return dict(tag=tag,nz=0,R=0,S=0,deficit=0,knobs=0,rank=0,failing=len(st.fails),score=st.score())
    R=set()
    for q in NZ: R|=atom_eqs[q]
    Rl=sorted(R)
    S=sorted(q for q in set().union(*[set(z for c,z in ev.eq_terms[e][2]) for e in Rl]) if atom_eqs[q]<=R)
    Sset=set(S); cands=set()
    for q in S:
        if q in fr.csup: cands.update(FR0[b] for b in _bits(fr.csup[q]))
    knobs=[]
    for Y in sorted(cands):
        if Y not in FREE: continue
        h=st.clone().set_free({Y:st.fv.get(Y,0)+1})
        d={q:h.av[q]-st.av[q] for q in h.av if h.av[q]!=st.av[q]}
        if d and all(q in Sset for q in d): knobs.append(Y)
    b=[inner(st,e) for e in Rl]; A=[]
    for Y in knobs:
        h=st.clone().set_free({Y:st.fv.get(Y,0)+1})
        A.append([inner(h,e)-b[i] for i,e in enumerate(Rl)])
    n=len(knobs)
    rk=rank_q([[A[j][i] for j in range(n)] for i in range(len(Rl))],n) if n else 0
    z0=sum(1 for x in b if x==0)
    outside=len([e for e in st.fails if e not in R])
    return dict(tag=tag,nz=len(NZ),R=len(Rl),S=len(S),deficit=len(Rl)-len(S),knobs=n,rank=rk,
                failing=len(Rl)-z0+outside,score=39033-(len(Rl)-z0+outside))
if __name__=='__main__':
    t0=time.time()
    cal=make(WIT); r=price(cal,'WITNESS {642,28730,29854,31864}')
    print('CALIBRATION: score=%d nz=%d |R|=%d |S|=%d deficit=%d knobs=%d rank=%d failing=%d  (%.1fs)'%(
        r['score'],r['nz'],r['R'],r['S'],r['deficit'],r['knobs'],r['rank'],r['failing'],time.time()-t0),flush=True)
    print('calibration reproduces 39,026:',r['score']==39026,flush=True)
