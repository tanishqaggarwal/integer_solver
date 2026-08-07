"""Final: accept a shift ONLY if it does not increase the global nonzero-atom count.
That guard subsumes every scoping question -- handle-carrying atoms, c==1 atoms, and atoms with
NO handle (which cannot absorb anything and must stay exactly zero) are all covered by it."""
import sys, pickle, collections, json, time, itertools, os
from math import gcd
src=open('/home/user/integer_solver/solve_lab/agentL_work/closeS3.py').read().split("def close(")[0]
exec(src)
def nzcount(vv):
    relift(vv); r=E.run(vv)
    return sum(1 for x in r if x)
def solve_group3(vv,V,w,gen,base):
    """candidates from the violated atoms' roots; accept only if global nonzero count improves."""
    Rs=[]
    for a in V:
        C=fitc(vv,a,w,gen)
        if C is None: return None
        rs=roots_c(C,abs(SL[a])//p)
        if not rs: return None
        Rs.append((rs,abs(SL[a])//p))
    for combo in itertools.product(*[r for r,_ in Rs]):
        t=crt_list([(r,c) for r,(_,c) in zip(combo,Rs)])
        if not t: continue
        old=vv[w]; vv[w]=old+p*t
        n=nzcount(vv)                       # GLOBAL guard, by direct recomputation
        if n<base: return t
        vv[w]=old
    return None
def close(S,tag,outer_max=14):
    v,isl,valn=assignment(set(S),ORIENT); v[24468]=T1; v[18956]=T2
    vv=[0]*NV
    for k,x in v.items(): vv[k]=x
    for rd in range(60):
        bad=relift(vv)
        if not bad: break
        r=E.run(vv); fx=0
        for a in bad:
            i=E.residx[a]; cur=r[i]; sm=abs(SL[a])
            if cur%p: continue
            imm=[q for q in vars_of(E.atoms[a]) if q in SHIFT]
            for w in imm+[q for q in atomvalvars[a] if q in SHIFT and q not in imm]:
                old=vv[w]; vv[w]=old+p; d=E.run(vv)[i]-cur; vv[w]=old
                if d==0: continue
                g=gcd(d,sm)
                if cur%g: continue
                mm=sm//g
                t=(-(cur//g))*pow((d//g)%mm,-1,mm)%mm if mm>1 else 0
                vv[w]=old+p*t; fx+=1; break
        if fx==0: break
    gen=0
    for outer in range(outer_max):
        base=nzcount(vv); r=E.run(vv); gen+=1
        viol=[a for a in SL if r[E.residx[a]]!=0 and SL[a] and r[E.residx[a]]%abs(SL[a])!=0]
        if not viol: break
        wires=collections.defaultdict(list)
        for a in viol:
            for w in (set(q for q in vars_of(E.atoms[a]) if q in SHIFT)|
                      set(q for q in atomvalvars[a] if q in SHIFT)): wires[w].append(a)
        prog=0
        for w,ats in sorted(wires.items(),key=lambda kv:-len(kv[1])):
            V=[a for a in ats if influences(vv,a,w)]
            if not V: continue
            t=solve_group3(vv,V,w,gen,base)
            if t: prog+=1; base=nzcount(vv); gen+=1
        if prog==0: break
    relift(vv); r=E.run(vv)
    nz=[E.res[i] for i,x in enumerate(r) if x]
    json.dump({'x_%d'%i:vv[i] for i in range(NV) if vv[i]},open('close_%s.json'%tag,'w'))
    return nz
if __name__=='__main__':
    tag=sys.argv[1]; n=int(sys.argv[2])
    import random
    rnd=random.Random(7)
    S=[24601,2081] if n==2 else rnd.sample(M['live'],n)
    t0=time.time(); nz=close(S,tag); el=time.time()-t0
    print('|S|=%-3d %-4s  NONZERO ATOMS = %d of 9032   WALL = %.1f s  -> close_%s.json'%(
        len(S),tag,len(nz),el,tag),flush=True)
    for a in nz[:8]: print('    ',a[:110],flush=True)
