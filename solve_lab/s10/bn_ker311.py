"""bn_ker311: exact nullspace of the 311-atom maximal support, then the cone
question re-posed in nullspace coordinates (d variables, 311 sign constraints)."""
import os, sys, json, time
from fractions import Fraction
HERE='/home/user/integer_solver/solve_lab/s10'
sys.path.insert(0, HERE)
import bn_lib as B
import lib as L, tools as T, ad
import bn_kernel as K

bools=B.bools_map()
S=sorted(json.load(open(os.path.join(HERE,'bn_keycomp.json')))['keycomp'])
S,E,rows = K.build(S)
print(f'block: {len(S)} atoms, {len(E)} eqs, {len(rows)} nonempty rows',flush=True)
t0=time.time()
rank, basis, freec = K.nullspace(S, rows)
print(f'exact rank {rank}, nullity {len(S)-rank}  ({time.time()-t0:.0f}s)',flush=True)
# sign-adjust: t_a = val_a / c_a  ->  column scaling
d=len(basis)
Bm=[[basis[j][i]/bools[S[i]][1] for j in range(d)] for i in range(len(S))]   # 311 x d
json.dump({'rank':rank,'nullity':d,'S':S}, open(os.path.join(HERE,'bn_ker311.json'),'w'))
print('nullity =',d,flush=True)

# LP in nullspace coords: find y in R^d with By >= 0, sum(By)=1
# Phase-1 with free variables y = y+ - y-
n=2*d
rowsL=[]
for i in range(len(S)):
    r=[Bm[i][j] for j in range(d)]+[-Bm[i][j] for j in range(d)]
    rowsL.append(r)                     # (By)_i >= 0  -> subtract slack
tot=[sum(rowsL[i][j] for i in range(len(S))) for j in range(n)]
print(f'reduced LP: {len(S)} sign constraints + 1 normalisation, {n} vars',flush=True)
# standard form: rowsL[i] . z - s_i = 0 ; tot . z = 1 ; z,s >= 0
m=len(S)+1
N=n+len(S)+m
Tb=[]
for i in range(len(S)):
    r=[Fraction(x) for x in rowsL[i]]+[Fraction(-1) if j==i else Fraction(0) for j in range(len(S))]
    Tb.append(r)
Tb.append([Fraction(x) for x in tot]+[Fraction(0)]*len(S))
b=[Fraction(0)]*len(S)+[Fraction(1)]
for i in range(m):
    Tb[i]=Tb[i]+[Fraction(1) if j==i else Fraction(0) for j in range(m)]+[b[i]]
basisv=[n+len(S)+i for i in range(m)]
cost=[Fraction(0)]*(N+1)
for i in range(m):
    for j in range(N+1): cost[j]+=Tb[i][j]
for i in range(m): cost[n+len(S)+i]-=1
it=0; t0=time.time(); bland=False
while True:
    it+=1
    if time.time()-t0>900: print('reduced LP timeout',it,flush=True); sys.exit()
    if it>3000 and not bland: bland=True
    piv=-1
    if bland:
        for j in range(N):
            if cost[j]>0: piv=j; break
    else:
        bv=0
        for j in range(N):
            if cost[j]>bv: bv=cost[j]; piv=j
    if piv<0: break
    ratio=None; pr=-1
    for i in range(m):
        if Tb[i][piv]>0:
            r=Tb[i][N]/Tb[i][piv]
            if ratio is None or r<ratio or (r==ratio and basisv[i]<basisv[pr]): ratio=r; pr=i
    if pr<0: print('UNBOUNDED'); break
    pv=Tb[pr][piv]; Tb[pr]=[x/pv for x in Tb[pr]]; rp=Tb[pr]
    for i in range(m):
        f=Tb[i][piv]
        if i!=pr and f: Tb[i]=[Tb[i][j]-f*rp[j] for j in range(N+1)]
    f=cost[piv]
    if f: cost=[cost[j]-f*rp[j] for j in range(N+1)]
    basisv[pr]=piv
obj=sum(Tb[i][N] for i in range(m) if basisv[i]>=n+len(S))
print(f'reduced LP: {it} pivots, artificial residual {obj} ({time.time()-t0:.0f}s)',flush=True)
print('RESULT:', 'TRIVIAL - no nonneg boolean kernel' if obj!=0 else 'FEASIBLE',flush=True)
json.dump({'rank':rank,'nullity':d,'status':('TRIVIAL' if obj!=0 else 'FEASIBLE')},
          open(os.path.join(HERE,'bn_ker311.json'),'w'))
