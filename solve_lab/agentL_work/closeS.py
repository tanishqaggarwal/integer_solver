"""Fixed solve_group (preserve already-satisfied c>1 atoms), NONZERO-ATOM metric, dump + check."""
import sys, pickle, collections, json, time, itertools
from math import gcd
src=open('/home/user/integer_solver/solve_lab/agentL_work/solve927g.py').read().split("#MAINSTART")[0]
exec(src)
CGTSET={a for a,s in SL.items() if s and s%p==0 and abs(s)//p>1}
def build_and_close(S,tag,outer_max=10):
    v,isl,valn=assignment(set(S),ORIENT); v[24468]=T1; v[18956]=T2
    vv=[0]*NV
    for k,x in v.items(): vv[k]=x
    for rd in range(60):                      # greedy first
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
    for outer in range(outer_max):
        r=E.run(vv)
        viol=[a for a in CGTSET if r[E.residx[a]]%(abs(SL[a]))!=0 and r[E.residx[a]]%p==0]
        if not viol: break
        cand=collections.defaultdict(list)
        for a in viol:
            ws=[q for q in vars_of(E.atoms[a]) if q in SHIFT]
            ws+=[q for q in atomvalvars[a] if q in SHIFT and q not in ws]
            for w in ws:
                if influences(vv,a,w): cand[w].append(a)
        prog=0
        for w,ats in sorted(cand.items(),key=lambda kv:-len(kv[1])):
            # THE S6i FIX: also require every ALREADY-SATISFIED c>1 atom that w touches to stay satisfied
            keep=[a for a in CGTSET if a not in ats and influences(vv,a,w)
                  and E.run(vv)[E.residx[a]]%abs(SL[a])==0]
            t=solve_group(vv,ats+keep,w)
            if t is not None and t!=0:
                vv[w]+=p*t; prog+=1
        if prog==0: break
    r=E.run(vv)
    nz=[E.res[i] for i,x in enumerate(r) if x]
    json.dump({'x_%d'%i:vv[i] for i in range(NV) if vv[i]},open('close_%s.json'%tag,'w'))
    print('|S|=%-3d  NONZERO ATOMS = %d of 9032   -> close_%s.json'%(len(S),len(nz),tag),flush=True)
    for a in nz[:6]: print('        ',a[:110],flush=True)
    return len(nz)
if __name__=='__main__':
    import random
    rnd=random.Random(7)
    todo=[(rnd.sample(M['live'],3),'S3'),(rnd.sample(M['live'],5),'S5'),
          (rnd.sample(M['live'],8),'S8'),([24601,2081],'S2'),(rnd.sample(M['live'],17),'S17')]
    for S,tag in todo:
        t0=time.time()
        try: build_and_close(S,tag)
        except Exception as e: print('|S|=%d %s FAILED %r'%(len(S),tag,e),flush=True)
        print('        (%.0fs)'%(time.time()-t0),flush=True)
