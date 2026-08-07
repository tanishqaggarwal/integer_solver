"""Stage B on the rank>deficit detach sets: EXACT integer optimum over the knob lattice."""
import ev, json, time, itertools, sys
from detach import make, POOL, WIT, fr, inner, _bits, rank_q, FREE, FR0, atom_eqs
from fractions import Fraction
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
def optimum(D,verbose=True):
    st=make(list(D))
    NZ=set(st.nz())
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
    # joint-linearity check
    g=st.clone()
    for Y in knobs: g.set_free({Y:st.fv.get(Y,0)+1})
    lin=all(inner(g,e)==b[i]+sum(A[j][i] for j in range(n)) for i,e in enumerate(Rl))
    rk=rank_q([[A[j][i] for j in range(n)] for i in range(len(Rl))],n) if n else 0
    z0=sum(1 for x in b if x==0)
    outside=len([e for e in st.fails if e not in R])
    best=z0; bt=None
    for size in range(min(n,rk,8),0,-1):
        if size<=best: break
        cnt=0
        for rows in itertools.combinations(range(len(Rl)),size):
            cnt+=1
            if cnt>400000: break
            t=solve_int(A,b,list(rows),n)
            if t is None: continue
            z=sum(1 for i in range(len(Rl)) if b[i]+sum(t[j]*A[j][i] for j in range(n))==0)
            if z>best: best=z; bt=t
    fail=len(Rl)-best+outside
    if verbose:
        print('D=%-22s |R|=%-3d |S|=%-3d deficit=%-3d knobs=%-3d rank=%-3d linear=%s z0=%-3d optimum=%-3d outside=%d failing=%d score=%d'%(
            str(list(D)),len(Rl),len(S),len(Rl)-len(S),n,rk,lin,z0,best,outside,fail,39033-fail),flush=True)
    return fail,bt,knobs,st,39033-fail
if __name__=='__main__':
    print('--- calibration ---',flush=True)
    optimum(WIT)
    print('--- rank>deficit candidates ---',flush=True)
    cands=[[17499]]+[[17499,v] for v in (20492,21279,22665,23642,26064,28599,28961,29854,31864,33347)]
    bestscore=39026; bestpack=None
    for D in cands:
        f,t,kn,st,sc=optimum(D)
        if sc>bestscore:
            bestscore=sc
            g=st.clone()
            if t:
                for j,Y in enumerate(kn):
                    if t[j]: g.set_free({Y:st.fv.get(Y,0)+t[j]})
            bestpack=(D,g,sc)
            print('  *** BEATS 39026: %d ***'%sc,flush=True)
    if bestpack:
        D,g,sc=bestpack
        out='H_%d_detach.json'%sc
        json.dump({('x_%d'%i):g.v[i] for i in range(38748) if g.v[i]!=0},open(out,'w'))
        print('WROTE',out,flush=True)
    else:
        print('\nno detach set beat 39,026 under the exact optimum.',flush=True)
