import sys, pickle, random, collections, re, json, time
src=open('/home/user/integer_solver/solve_lab/agentV_work/mirror/L/calib2.py').read().split("# numeric perm repair")[0]
exec(src)
defrhs={c[1]:c[2] for c in (E.cls[a] for a in E.order)}
C2=pickle.load(open('calib2.pkl','rb')); perm.update(C2['perm']); ORIENT=C2['ORIENT']
T1=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
T2=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
# handle var per atom (free var whose variation moves the atom by a multiple of p) -> h = V/(M*p)
Hd=pickle.load(open('handles.pkl','rb')); handle=set(Hd['handle'])
from circ2 import vars_of
freeall={}
def fa(v):
    if v in freeall: return freeall[v]
    if v not in defrhs: freeall[v]={v}; return freeall[v]
    freeall[v]=set(); s=set()
    for u in vars_of(defrhs[v]): s|=fa(u)
    freeall[v]=s; return s
sys.setrecursionlimit(100000)
atomh={}
for a in E.res:
    s=set()
    for u in vars_of(E.atoms[a]): s|=fa(u)
    hs=[u for u in s if u in handle]
    atomh[a]=hs
print('atoms with exactly 1 handle:',sum(1 for a in E.res if len(atomh[a])==1),'of',len(E.res))

def full(S,extra=None):
    v,isl,valn=assignment(set(S),ORIENT)
    v[24468]=T1; v[18956]=T2
    if extra: v.update(extra)
    vv,r=run(v)
    return vv,r,valn

def lift(vv):
    """set every handle var so its atom vanishes exactly over Z; returns #unliftable"""
    bad=[]
    for a in E.res:
        hs=atomh[a]
        if len(hs)!=1: continue
        h=hs[0]
        i=E.residx[a]
        v0=vv[h]; vv[h]=v0
        # atom is affine in h with slope s
        cur=E.run(vv)[i]
        vv[h]=v0+1; s=E.run(vv)[i]-cur; vv[h]=v0
        if s==0:
            if cur: bad.append((a,'noslope'))
            continue
        if cur % s: bad.append((a,'nondiv',cur,s)); continue
        vv[h]=v0-cur//s
    return bad

if __name__=='__main__':
    t0=time.time()
    S=[M['live'][0]]
    vv,r,valn=full(S)
    nz=[i for i,x in enumerate(r) if x%p]
    print('|S|=1  nonzero atoms mod p after target pins:',len(nz))
    for i in nz: print('   ',E.res[i][:150])
    bad=E.score(r); print('failing equations (F-model score):',len(bad),'=> ',39033-len(bad))
    print('root fold',valn[ROOT],' T1modp',T1%p,' T2modp',T2%p)
