"""Sweep |S|=2,3,5,8,17 with the S6i fix; influence map built STRUCTURALLY (no probing)."""
import sys, pickle, collections, json, time, itertools
from math import gcd
src=open('/home/user/integer_solver/solve_lab/agentL_work/solve927g.py').read().split("#MAINSTART")[0]
exec(src)
CGT2={a for a,s in SL.items() if s and s%p==0 and abs(s)//p>1}
# ---- structural influence map, built ONCE ----
W2A=collections.defaultdict(list)
for a in CGT2:
    ws=set(q for q in vars_of(E.atoms[a]) if q in SHIFT)
    ws|=set(q for q in atomvalvars[a] if q in SHIFT)
    for w in ws: W2A[w].append(a)
print('structural influence map: %d wires, mean %.1f c>1 atoms/wire, max %d'%(
    len(W2A), sum(len(v) for v in W2A.values())/max(1,len(W2A)), max(len(v) for v in W2A.values())),flush=True)
def close(S,tag,outer_max=12):
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
    for outer in range(outer_max):
        r=E.run(vv)
        viol=[a for a in CGT2 if r[E.residx[a]]!=0 and r[E.residx[a]]%abs(SL[a])!=0]
        if not viol: break
        wires=collections.defaultdict(list)
        for a in viol:
            for w in (set(q for q in vars_of(E.atoms[a]) if q in SHIFT)|
                      set(q for q in atomvalvars[a] if q in SHIFT)): wires[w].append(a)
        prog=0
        for w,ats in sorted(wires.items(),key=lambda kv:-len(kv[1])):
            ats=[a for a in ats if influences(vv,a,w)]
            if not ats: continue
            keep=[a for a in W2A[w] if a not in ats and influences(vv,a,w)
                  and r[E.residx[a]]%abs(SL[a])==0]
            t=solve_group(vv,ats+keep,w)
            if t: vv[w]+=p*t; prog+=1; r=E.run(vv)
        if prog==0: break
    r=E.run(vv)
    nz=[E.res[i] for i,x in enumerate(r) if x]
    json.dump({'x_%d'%i:vv[i] for i in range(NV) if vv[i]},open('close_%s.json'%tag,'w'))
    print('|S|=%-3d %-4s NONZERO ATOMS = %d of 9032  -> close_%s.json'%(len(S),tag,len(nz),tag),flush=True)
    for a in nz[:5]: print('       ',a[:110],flush=True)
    return len(nz)
if __name__=='__main__':
    import random
    rnd=random.Random(7)
    for S,tag in [([24601,2081],'S2'),(rnd.sample(M['live'],3),'S3'),(rnd.sample(M['live'],5),'S5'),
                  (rnd.sample(M['live'],8),'S8'),(rnd.sample(M['live'],17),'S17')]:
        t0=time.time()
        try: close(S,tag)
        except Exception as e: print('|S|=%d %s FAILED %r'%(len(S),tag,e),flush=True)
        print('       (%.0fs)'%(time.time()-t0),flush=True)
