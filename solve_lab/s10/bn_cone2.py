"""bn_cone2: faster exact cone decision (Dantzig pricing + Bland fallback).

Decides whether {t >= 0, M t = 0, sum t = 1} is feasible, where
M[e][a] = coeff(e,a) * c_a  and t_a = x_a(x_a - 1) >= 0.

usage: bn_cone2.py <which>   which in {free, blk29, allcore}
Writes result to bn_cone2_<which>.json
"""
import os, sys, json, collections, time
from fractions import Fraction
HERE='/home/user/integer_solver/solve_lab/s10'
sys.path.insert(0, HERE)
import bn_lib as B
import lib as L, tools as T, ad

bools=B.bools_map()

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
        if verbose: print(f'   joint-peel r{it}: -{len(kill)} -> {len(cand)}',flush=True)
    return cand

def cone(S, tlimit=1500):
    S=sorted(S); n=len(S)
    if n==0: return 'EMPTY', None
    idx={a:i for i,a in enumerate(S)}
    E=sorted(set().union(*[set(L.atom2eq[a]) for a in S]))
    rows=[]
    for e in E:
        m,sq,co=L.eq_atoms[e]
        r=[0]*n; any_=False
        for a in co:
            if a in idx: r[idx[a]]=co[a]*bools[a][1]; any_=True
        if any_: rows.append(r)
    rows.append([1]*n)
    b=[Fraction(0)]*(len(rows)-1)+[Fraction(1)]
    m=len(rows); N=n+m
    print(f'   LP {m} x {n}',flush=True)
    Tb=[[Fraction(x) for x in rows[i]]+[Fraction(1) if j==i else Fraction(0) for j in range(m)]+[b[i]]
        for i in range(m)]
    basis=[n+i for i in range(m)]
    cost=[Fraction(0)]*(N+1)
    for i in range(m):
        for j in range(N+1): cost[j]+=Tb[i][j]
    for i in range(m): cost[n+i]-=1
    it=0; t0=time.time(); bland=False
    while True:
        it+=1
        if time.time()-t0>tlimit: return 'TIMEOUT', it
        if it>4000 and not bland: bland=True
        piv=-1
        if bland:
            for j in range(N):
                if cost[j]>0: piv=j; break
        else:
            bestv=0
            for j in range(N):
                if cost[j]>bestv: bestv=cost[j]; piv=j
        if piv<0: break
        ratio=None; pr=-1
        for i in range(m):
            if Tb[i][piv]>0:
                r=Tb[i][N]/Tb[i][piv]
                if ratio is None or r<ratio or (r==ratio and basis[i]<basis[pr]):
                    ratio=r; pr=i
        if pr<0: return 'UNBOUNDED', it
        pv=Tb[pr][piv]
        Tb[pr]=[x/pv for x in Tb[pr]]
        rp=Tb[pr]
        for i in range(m):
            f=Tb[i][piv]
            if i!=pr and f:
                Ti=Tb[i]
                Tb[i]=[Ti[j]-f*rp[j] for j in range(N+1)]
        f=cost[piv]
        if f: cost=[cost[j]-f*rp[j] for j in range(N+1)]
        basis[pr]=piv
    obj=sum(Tb[i][N] for i in range(m) if basis[i]>=n)
    print(f'   pivots {it}  artificial residual {obj}  ({time.time()-t0:.0f}s)',flush=True)
    if obj!=0: return 'TRIVIAL', None
    t={S[basis[i]]: Tb[i][N] for i in range(m) if basis[i]<n and Tb[i][N]}
    return 'FEASIBLE', t

if __name__=='__main__':
    which=sys.argv[1]
    src={'free':'bn_fcore.json','allcore':'bn_core.json'}
    if which=='blk29':
        cand=set(json.load(open(os.path.join(HERE,'bn_defic.json')))['all']['S'])
    else:
        cand=set(json.load(open(os.path.join(HERE,src[which]))))
    print(f'{which}: start {len(cand)}',flush=True)
    c=joint_peel(cand)
    print(f'  joint-peel survivors {len(c)}',flush=True)
    st,res=cone(c)
    print(f'  RESULT {which}: {st}',flush=True)
    out={'which':which,'peel':sorted(c),'status':st,
         'sol':({str(a):[v.numerator,v.denominator] for a,v in res.items()}
                if isinstance(res,dict) else None)}
    json.dump(out, open(os.path.join(HERE,f'bn_cone2_{which}.json'),'w'))
    print('saved',flush=True)
