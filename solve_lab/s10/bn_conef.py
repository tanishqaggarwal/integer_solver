"""bn_conef: fast float Phase-1 simplex for the cone question, on the 311-atom
maximal support.  Any FEASIBLE answer is then re-verified exactly on its support.
"""
import os, sys, json, collections, time
from fractions import Fraction
HERE='/home/user/integer_solver/solve_lab/s10'
sys.path.insert(0, HERE)
import bn_lib as B
import lib as L, tools as T, ad

bools=B.bools_map()
S=sorted(json.load(open(os.path.join(HERE,'bn_keycomp.json')))['keycomp'])
n=len(S); idx={a:i for i,a in enumerate(S)}
E=sorted(set().union(*[set(L.atom2eq[a]) for a in S]))
rows=[]
for e in E:
    m,sq,co=L.eq_atoms[e]
    r=[0.0]*n; any_=False
    for a in co:
        if a in idx: r[idx[a]]=float(co[a]*bools[a][1]); any_=True
    if any_: rows.append(r)
rows.append([1.0]*n)
b=[0.0]*(len(rows)-1)+[1.0]
m=len(rows); N=n+m
print(f'float LP {m} x {n}',flush=True)
Tb=[rows[i]+[1.0 if j==i else 0.0 for j in range(m)]+[b[i]] for i in range(m)]
basis=[n+i for i in range(m)]
cost=[0.0]*(N+1)
for i in range(m):
    for j in range(N+1): cost[j]+=Tb[i][j]
for i in range(m): cost[n+i]-=1.0
EPS=1e-9
t0=time.time(); it=0
while True:
    it+=1
    if time.time()-t0>1200: print('float LP timeout at',it); break
    piv=-1; bv=EPS
    for j in range(N):
        if cost[j]>bv: bv=cost[j]; piv=j
    if piv<0: break
    ratio=None; pr=-1
    for i in range(m):
        if Tb[i][piv]>EPS:
            r=Tb[i][N]/Tb[i][piv]
            if ratio is None or r<ratio-1e-12: ratio=r; pr=i
    if pr<0: print('unbounded'); break
    pv=Tb[pr][piv]; rp=[x/pv for x in Tb[pr]]; Tb[pr]=rp
    for i in range(m):
        f=Tb[i][piv]
        if i!=pr and abs(f)>1e-14:
            Ti=Tb[i]
            Tb[i]=[Ti[j]-f*rp[j] for j in range(N+1)]
    f=cost[piv]
    if abs(f)>1e-14: cost=[cost[j]-f*rp[j] for j in range(N+1)]
    basis[pr]=piv
obj=sum(Tb[i][N] for i in range(m) if basis[i]>=n)
print(f'pivots {it}  artificial residual {obj:.3e}  ({time.time()-t0:.0f}s)',flush=True)
if obj>1e-7:
    print('RESULT: cone is TRIVIAL (float) -> no nonneg boolean kernel on the 311 core')
    res={'status':'TRIVIAL_float','obj':obj}
else:
    sup={S[basis[i]]:Tb[i][N] for i in range(m) if basis[i]<n and Tb[i][N]>1e-9}
    print('RESULT: FEASIBLE (float), support',len(sup))
    print('  ',{f'a{a}':round(v,6) for a,v in list(sup.items())[:15]})
    res={'status':'FEASIBLE_float','support':{str(a):v for a,v in sup.items()}}
json.dump(res, open(os.path.join(HERE,'bn_conef.json'),'w'))
print('saved bn_conef.json',flush=True)
