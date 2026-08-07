"""Fixed solve_group: enumerate the VIOLATED atom's roots only; preserve `keep` atoms by
forcing t == 0 mod c_keep (t=0 is always a root of a currently-satisfied atom).  No enumeration
of keep root sets at all.  fit() cached per (atom,wire) within a pass.  Verified by recomputation."""
import sys, pickle, collections, json, time, itertools
from math import gcd
src=open('/home/user/integer_solver/solve_lab/agentT_work/mirror/L/solve927g.py').read().split("#MAINSTART")[0]
exec(src)
CGT2={a for a,s in SL.items() if s and s%p==0 and abs(s)//p>1}
W2A=collections.defaultdict(list)
for a in CGT2:
    ws=set(q for q in vars_of(E.atoms[a]) if q in SHIFT)|set(q for q in atomvalvars[a] if q in SHIFT)
    for w in ws: W2A[w].append(a)
FITC={}
def fitc(vv,a,w,gen):
    k=(a,w,gen)
    if k not in FITC: FITC[k]=fit(vv,E.residx[a],w)[0]
    return FITC[k]
def lcm(x,y): return x//gcd(x,y)*y
def roots_c(C,c):
    """all residues t mod c with poly(t)==0 mod c, via prime powers + CRT"""
    per=[]
    for q,e in factor(c).items():
        rs=sorted(rootset_pp(C,q,e))
        if not rs: return []
        per.append((rs,q**e))
    out=[]
    for combo in itertools.product(*[r for r,_ in per]):
        t=crt_list([(r,m) for r,(_,m) in zip(combo,per)])
        if t is not None: out.append(t)
    return out
def solve_group2(vv,V,K,w,gen):
    Rs=[]
    for a in V:
        C=fitc(vv,a,w,gen)
        if C is None: return None
        c=abs(SL[a])//p
        rs=roots_c(C,c)
        if not rs: return None
        Rs.append((rs,c))
    Lk=1
    for a in K: Lk=lcm(Lk,abs(SL[a])//p)
    for combo in itertools.product(*[r for r,_ in Rs]):
        pairs=[(r,c) for r,(_,c) in zip(combo,Rs)]
        if K: pairs.append((0,Lk))
        t=crt_list(pairs)
        if t is None or t==0: continue
        ok=True
        for a in V+K:                                  # GUARD: direct recomputation
            if probe(vv,E.residx[a],[w],[t])%abs(SL[a])!=0: ok=False; break
        if ok: return t
    return None
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
    gen=0
    for outer in range(outer_max):
        r=E.run(vv); gen+=1
        viol=[a for a in CGT2 if r[E.residx[a]]!=0 and r[E.residx[a]]%abs(SL[a])!=0]
        if not viol: break
        wires=collections.defaultdict(list)
        for a in viol:
            for w in (set(q for q in vars_of(E.atoms[a]) if q in SHIFT)|
                      set(q for q in atomvalvars[a] if q in SHIFT)): wires[w].append(a)
        prog=0
        for w,ats in sorted(wires.items(),key=lambda kv:-len(kv[1])):
            V=[a for a in ats if influences(vv,a,w)]
            if not V: continue
            K=[a for a in W2A[w] if a not in V and r[E.residx[a]]%abs(SL[a])==0 and influences(vv,a,w)]
            t=solve_group2(vv,V,K,w,gen)
            if t: vv[w]+=p*t; prog+=1; r=E.run(vv); gen+=1
        if prog==0: break
    r=E.run(vv)
    nz=[E.res[i] for i,x in enumerate(r) if x]
    json.dump({'x_%d'%i:vv[i] for i in range(NV) if vv[i]},open('close_%s.json'%tag,'w'))
    return nz
if __name__=='__main__':
    tag=sys.argv[1]; n=int(sys.argv[2])
    import random
    rnd=random.Random(7)
    S=[24601,2081] if n==2 else rnd.sample(M['live'],n)
    t0=time.time(); nz=close(S,tag); el=time.time()-t0
    print('|S|=%-3d %-4s  NONZERO ATOMS = %d of 9032   WALL CLOCK = %.1f s  -> close_%s.json'%(
        len(S),tag,len(nz),el,tag),flush=True)
    for a in nz[:6]: print('    ',a[:110],flush=True)
