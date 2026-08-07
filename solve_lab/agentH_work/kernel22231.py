"""Rank-raising sweep, stage 2: is there an integer DIRECTION that moves a22231 while every atom
   outside the 8-atom region stays put?  That is the only remaining way to get +0 new equations."""
import frameB as FB, ev, json, time, math
from frameB import Frame, State
from collections import defaultdict
from fractions import Fraction
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
atom_eqs=defaultdict(set)
for i,(m,sq,tl) in enumerate(ev.eq_terms):
    for c,a in tl: atom_eqs[a].add(i)
fr=Frame([642,28730,29854,31864])
W=json.load(open('../best/new_instance_partial_39026.json'))
v=[0]*38748
for k,val in W.items(): v[int(k[2:]) if k.startswith('x_') else int(k)]=int(val)
base=State(fr,{u:v[u] for u in fr.free if v[u]!=0})
IN8={22229,22230,22231,35758,35759,35760,35761,35762}
MOVERS=[2081,4287,4432,12553,28730]
# widen: every free input that moves any atom the movers disturb outside IN8
coll=set()
for X in MOVERS:
    g=base.clone().set_free({X:base.fv.get(X,0)+1})
    coll |= set(a for a in g.av if g.av[a]!=base.av[a] and a not in IN8)
print('collateral atoms to cancel:',sorted(coll))
helpers=set(MOVERS)
for a in coll: helpers |= set(fr.SUPV.get(a,[]))
helpers={X for X in helpers if X in set(fr.free)}
print('candidate knobs (movers + anything touching the collateral):',len(helpers))
# measure the derivative of every atom (outside IN8, plus a22231) wrt each candidate knob
def deltas(X,d=1):
    g=base.clone().set_free({X:base.fv.get(X,0)+d})
    return {a:g.av[a]-base.av[a] for a in g.av if g.av[a]!=base.av[a]}
rowsq={}   # squares handled by their linear root
D={}
allatoms=set()
nonlin=[]
for X in sorted(helpers):
    d1=deltas(X,1); d2=deltas(X,2)
    lin=all(d2.get(a,0)==2*d1[a] for a in d1) and set(d2)==set(d1)
    D[X]=(d1,lin)
    if not lin: nonlin.append(X)
    allatoms|=set(d1)
print('candidate knobs with NONLINEAR atom response:',len(nonlin),nonlin[:10])
OUT=sorted(a for a in allatoms if a not in IN8)
print('atoms outside the region that must be held fixed:',len(OUT))
# build the rational matrix: rows = OUT atoms (using linear roots for squares), cols = knobs
KN=sorted(D)
def entry(X,a):
    d1,lin=D[X]
    val=d1.get(a,0)
    if lin: return Fraction(val)
    # square atom: value is (linear)^2 ; use the signed root
    if val==0: return Fraction(0)
    r=math.isqrt(abs(val))
    if r*r==abs(val): return Fraction(r if val>0 else -r)
    return None
M=[]; skipped=[]
for a in OUT:
    row=[entry(X,a) for X in KN]
    if any(x is None for x in row): skipped.append(a); continue
    M.append((a,row))
print('rows usable %d, skipped (no exact root) %d %s'%(len(M),len(skipped),skipped[:6]))
targ=[Fraction(D[X][0].get(22231,0)) for X in KN]
print('knobs with nonzero a22231 derivative:',[KN[j] for j in range(len(KN)) if targ[j]!=0])
# solve: find t with M t = 0 and targ . t != 0   <=>  targ NOT in rowspace(M)
A=[r for _,r in M]+[targ]
def rank(rows,ncols):
    rr=[list(r) for r in rows]; piv=0
    for c in range(ncols):
        k=None
        for i in range(piv,len(rr)):
            if rr[i][c]!=0: k=i; break
        if k is None: continue
        rr[piv],rr[k]=rr[k],rr[piv]
        pv=rr[piv][c]
        for i in range(len(rr)):
            if i!=piv and rr[i][c]!=0:
                f=rr[i][c]/pv
                rr[i]=[rr[i][j]-f*rr[piv][j] for j in range(ncols)]
        piv+=1
        if piv==len(rr): break
    return piv
n=len(KN)
rM=rank([r for _,r in M],n); rA=rank(A,n)
print('rank(collateral rows) = %d   rank(collateral rows + a22231 row) = %d   knobs = %d'%(rM,rA,n))
if rA>rM:
    print('=> a22231 is INDEPENDENT of the collateral rows: a zero-collateral direction EXISTS.')
else:
    print('=> a22231 lies IN the row space of the collateral constraints:')
    print('   every direction that moves a22231 necessarily moves something outside the region.')
    print('   NO knob and NO combination can achieve +0 new equations.  7 is the floor.')
