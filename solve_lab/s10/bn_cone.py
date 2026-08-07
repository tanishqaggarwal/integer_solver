"""bn_cone: does a NONZERO NONNEGATIVE kernel vector exist?

t_a = x(x-1) is ALWAYS >= 0 (image {0,2,6,12,20,...}).  So a realisable
configuration needs  t >= 0, t != 0, sum_a coeff(e,a)*c_a*t_a = 0 for all e.
Exact rational Phase-1 simplex decides this.

Peels applied first (both are sound implications, they never discard a solution):
  zero-peel: equation with exactly one surviving boolean atom => that t = 0
  sign-peel: equation whose surviving coeff*c_a all share a sign => all those t = 0
"""
import os, sys, json, collections
from fractions import Fraction
HERE='/home/user/integer_solver/solve_lab/s10'
sys.path.insert(0, HERE)
import bn_lib as B
import lib as L, tools as T, ad

bools = B.bools_map()
BA = set(bools)

def joint_peel(cand, verbose=True):
    cand=set(cand); it=0
    while True:
        it+=1
        eqs=collections.defaultdict(list)
        for a in cand:
            for e,co in L.atom2eq[a].items(): eqs[e].append((a, co*bools[a][1]))
        kill=set()
        for e,As in eqs.items():
            if len(As)==1: kill.add(As[0][0]); continue
            if len(set(1 if c>0 else -1 for _,c in As))==1:
                kill.update(a for a,_ in As)
        if not kill: break
        cand-=kill
        if verbose: print(f'   joint-peel r{it}: -{len(kill)} -> {len(cand)}', flush=True)
    return cand

# ---------------- exact rational Phase-1 simplex ----------------
def cone_nonzero(S, verbose=True):
    """Decide if {t>=0, Mt=0, sum t = 1} is feasible.  Returns t (dict) or None."""
    S=sorted(S); n=len(S)
    if n==0: return None
    idx={a:i for i,a in enumerate(S)}
    E=sorted(set().union(*[set(L.atom2eq[a]) for a in S]))
    rows=[]
    for e in E:
        m,sq,co=L.eq_atoms[e]
        r=[0]*n
        any_=False
        for a in co:
            if a in idx: r[idx[a]] = co[a]*bools[a][1]; any_=True
        if any_: rows.append(r)
    rows.append([1]*n)                       # sum t = 1
    b=[Fraction(0)]*(len(rows)-1)+[Fraction(1)]
    m=len(rows)
    if verbose: print(f'   LP: {m} rows x {n} cols', flush=True)
    # tableau: [A | I_art], basis = artificials
    T_=[[Fraction(x) for x in rows[i]]+[Fraction(1) if j==i else Fraction(0)
        for j in range(m)]+[b[i]] for i in range(m)]
    basis=[n+i for i in range(m)]
    N=n+m
    # cost row: minimize sum artificials
    cost=[Fraction(0)]*(N+1)
    for i in range(m):
        for j in range(N+1): cost[j]-=T_[i][j]
    for i in range(m): cost[n+i]+=Fraction(1)*0   # artificials have cost 1, reduced=0
    # recompute properly: z_j - c_j  with c_j = 1 for artificials, 0 else
    cost=[Fraction(0)]*(N+1)
    for i in range(m):
        for j in range(N+1): cost[j]+=T_[i][j]
    for i in range(m): cost[n+i]-=Fraction(1)
    it=0
    while True:
        it+=1
        if it>20000: print('   LP iteration cap'); return None
        piv=-1
        for j in range(N):                    # Bland's rule
            if cost[j]>0: piv=j; break
        if piv<0: break
        ratio=None; pr=-1
        for i in range(m):
            if T_[i][piv]>0:
                r=T_[i][N]/T_[i][piv]
                if ratio is None or r<ratio or (r==ratio and basis[i]<basis[pr]):
                    ratio=r; pr=i
        if pr<0: print('   LP unbounded'); return None
        pv=T_[pr][piv]
        T_[pr]=[x/pv for x in T_[pr]]
        for i in range(m):
            if i!=pr and T_[i][piv]:
                f=T_[i][piv]
                T_[i]=[T_[i][j]-f*T_[pr][j] for j in range(N+1)]
        if cost[piv]:
            f=cost[piv]
            cost=[cost[j]-f*T_[pr][j] for j in range(N+1)]
        basis[pr]=piv
    obj=sum(T_[i][N] for i in range(m) if basis[i]>=n)
    if verbose: print(f'   LP done in {it} pivots, artificial residual {obj}', flush=True)
    if obj!=0: return None
    t={}
    for i in range(m):
        if basis[i]<n and T_[i][N]: t[S[basis[i]]]=T_[i][N]
    return t

def istri(t):
    if t<0 or t%2: return False
    mm=1+4*t; r=int(mm**0.5)
    while r*r>mm: r-=1
    while (r+1)*(r+1)<=mm: r+=1
    return r*r==mm

def kfor(t):
    mm=1+4*t; r=int(mm**0.5)
    while r*r>mm: r-=1
    while (r+1)*(r+1)<=mm: r+=1
    return (1+r)//2 if r*r==mm else None

if __name__=='__main__':
    OUT=os.path.join(HERE,'bn_cone.json')
    res={}
    fcore=set(json.load(open(os.path.join(HERE,'bn_fcore.json'))))
    core=set(json.load(open(os.path.join(HERE,'bn_core.json'))))
    d=json.load(open(os.path.join(HERE,'bn_defic.json')))
    blk=set(d['all']['S'])
    for name,cand in (('free',fcore),('blk29',blk),('allcore',core)):
        print(f'--- {name}: start {len(cand)} ---',flush=True)
        c=joint_peel(cand)
        print(f'   joint-peel survivors {len(c)}',flush=True)
        res[name+'_peel']=sorted(c)
        json.dump(res,open(OUT,'w'))
        if not c:
            print(f'   {name}: EMPTY -> no nonzero boolean configuration possible')
            res[name+'_cone']=None; json.dump(res,open(OUT,'w')); continue
        t=cone_nonzero(c)
        if t is None:
            print(f'   {name}: cone is TRIVIAL -> every boolean atom forced to zero')
            res[name+'_cone']=None
        else:
            print(f'   {name}: FEASIBLE nonneg kernel, support {len(t)}')
            print('      ',{f'a{a}':str(v) for a,v in list(t.items())[:12]})
            res[name+'_cone']={str(a):[v.numerator,v.denominator] for a,v in t.items()}
        json.dump(res,open(OUT,'w'))
    print('saved',OUT)
