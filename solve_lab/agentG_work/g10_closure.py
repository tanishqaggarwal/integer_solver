"""Close the non-boolean symbol set: symbolize, find symbolic checks, take their
free-input supports, add new non-boolean ones, repeat until fixed point."""
import os, sys, json, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym
from gsym import *
import suppfree
src = sys.argv[1] if len(sys.argv)>1 else '/home/user/integer_solver/solve_lab/s10/AG_39013.json'
v = L.load(src); ad.fwd(v, rounds=6); vm=[x%P for x in v]
idx, freelist, vs = suppfree.build(vm, modp=True)
FREESET=set(freelist)
def isbool(u):
    for a in L.var_atoms[u]:
        pl=L.polys[a]
        if len(pl)==2:
            ks=list(pl.keys())
            if sorted(map(len,ks))==[1,2] and (u,) in pl and (u,u) in pl and pl[(u,)]==-pl[(u,u)]:
                return True
    return False
BOOLCACHE={}
def B(u):
    if u not in BOOLCACHE: BOOLCACHE[u]=isbool(u)
    return BOOLCACHE[u]
def supp(a):
    m=0
    for w in L.avars[a]: m |= vs[w] if w<len(vs) else 0
    s={freelist[i] for i in range(len(freelist)) if (m>>i)&1}
    s |= {w for w in L.avars[a] if w in FREESET}
    return s
S=set(json.load(open('supp8.json'))['nonbool'])
for it in range(12):
    SY=sorted(S)
    val=gsym.build(v,SY,cap=None,verbose=False)
    n=len(SY)
    symchecks=[]
    for a in gsym.check_atoms():
        f=gsym.evalpoly_sym(a,val,n,None)
        if not isinstance(f,int): symchecks.append(a)
        elif f%P: symchecks.append(a)   # nonzero constant - also a target
    new=set()
    for a in symchecks:
        for u in supp(a):
            if not B(u): new.add(u)
    print('iter %d: |S|=%d symchecks=%d newnonbool=%d' % (it,len(S),len(symchecks),len(new-S)))
    if new <= S: break
    S |= new
print('CLOSED non-boolean symbol set (%d): %s' % (len(S),sorted(S)))
json.dump(sorted(S), open('closed_nonbool.json','w'))
