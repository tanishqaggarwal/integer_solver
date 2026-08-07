"""solve_group: clear ALL conditions carried by a contended wire at once.
Root sets are intersected prime-power by prime-power, then CRT'd. Verified by recomputation."""
import sys, pickle, collections, json, time, itertools
from math import gcd
src=open('/home/user/integer_solver/solve_lab/agentT_work/mirror/L/solve927.py').read().split("if __name__")[0]
exec(src)
def rootset_pp(coeffs,q,e):
    m=q**e
    return {t for t in range(m) if peval(coeffs,t,m)==0}
def crt_list(pairs):
    t,mod=0,1
    for r,m in pairs:
        g=gcd(mod,m)
        if (r-t)%g: return None
        l=mod//g*m
        t=(t+mod*(((r-t)//g)*pow(mod//g,-1,m//g)%(m//g)))%l
        mod=l
    return t
def solve_group(vv,atoms,wire):
    """find t for `wire` clearing every atom in `atoms`; None if impossible."""
    per={}                       # prime -> (maxexp, set of allowed residues mod q^maxexp)
    for a in atoms:
        i=E.residx[a]; c=abs(SL[a])//p
        coeffs,td=fit(vv,i,wire)
        if coeffs is None: return None
        for q,e in factor(c).items():
            rs=rootset_pp(coeffs,q,e)
            if not rs: return None
            if q in per:
                e0,s0=per[q]
                if e>=e0:
                    s0={r for r in range(q**e) if r%(q**e0) in s0}
                    per[q]=(e,s0&rs)
                else:
                    rs={r for r in range(q**e0) if r%(q**e) in rs}
                    per[q]=(e0,s0&rs)
            else: per[q]=(e,rs)
            if not per[q][1]: return None
    keys=sorted(per)
    for combo in itertools.product(*[sorted(per[k][1]) for k in keys]):
        t=crt_list([(r,k**per[k][0]) for r,k in zip(combo,keys)])
        if t is None: continue
        ok=True                                  # GUARD: direct recomputation
        for a in atoms:
            i=E.residx[a]; c=abs(SL[a])//p
            if probe(vv,i,[wire],[t])%(c*p)!=0: ok=False; break
        if ok: return t
    return None
def influences(vv,a,w):
    i=E.residx[a]
    return probe(vv,i,[w],[1])!=probe(vv,i,[w],[0])
#MAINSTART
if __name__=='__main__':
    import random
    rnd=random.Random(7)
    S=rnd.sample(M['live'],17)
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
    print('=== |S|=17 JOINT solve ===',flush=True)
    for outer in range(15):
        left=relift(vv); r=E.run(vv)
        stuck=[a for a in left if r[E.residx[a]]%p==0]
        if not stuck:
            print('ALL DISCHARGED after %d outer rounds'%outer,flush=True); break
        print('round %d: %d stuck'%(outer,len(stuck)),flush=True)
        # candidate wires, and which stuck atoms each influences
        cand=collections.defaultdict(list)
        for a in stuck:
            ws=[q for q in vars_of(E.atoms[a]) if q in SHIFT]
            ws+=[q for q in atomvalvars[a] if q in SHIFT and q not in ws]
            for w in ws:
                if influences(vv,a,w): cand[w].append(a)
        done=set(); prog=0
        for w,ats in sorted(cand.items(),key=lambda kv:-len(kv[1])):
            ats=[a for a in ats if a not in done]
            if not ats: continue
            t0=time.time(); t=solve_group(vv,ats,w)
            if t is not None:
                vv[w]+=p*t; done.update(ats); prog+=1
                print('   wire x%-6d cleared %d condition(s) jointly, t=%-10d (%.0fs, verified)'%(
                    w,len(ats),t,time.time()-t0),flush=True)
        if prog==0:
            print('   no wire could clear its group -> RESIDUE, not oscillation',flush=True); break
    left=relift(vv); r=E.run(vv)
    stuck=[a for a in left if r[E.residx[a]]%p==0]
    nz=sum(1 for x in r if x)
    print('FINAL |S|=17: %d undischarged, %d nonzero atoms'%(len(stuck),nz),flush=True)
    for a in stuck: print('   residue c=%-10d %s'%(abs(SL[a])//p,a[:70]),flush=True)
