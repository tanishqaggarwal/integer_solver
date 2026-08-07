"""Bivariate solve for the cross-wire residue + component sizes of the shares-a-condition graph."""
import sys, pickle, collections, json, time, itertools, os
from math import gcd
src=open('/home/user/integer_solver/solve_lab/agentL_work/solve927g.py').read().split("if __name__")[0]
exec(src)
MYPID=os.getpid()
def fit2(vv,i,w1,w2,deg=4):
    """bivariate Newton fit of R(t1,t2)/p on a (deg+1)^2 grid"""
    G=[]
    for a in range(deg+1):
        row=[]
        for b in range(deg+1):
            y=probe(vv,i,[w1,w2],[a,b])
            if y%p: return None
            row.append(y//p)
        G.append(row)
    # forward differences in b, then in a
    D=[[r[:] for r in G]]
    for r in range(len(G)):
        d=[G[r][:]]
        for k in range(deg): d.append([d[k][j+1]-d[k][j] for j in range(len(d[k])-1)])
        G[r]=[d[k][0] for k in range(deg+1)]
    C=[]
    for c in range(deg+1):
        col=[G[r][c] for r in range(deg+1)]
        d=[col[:]]
        for k in range(deg): d.append([d[k][j+1]-d[k][j] for j in range(len(d[k])-1)])
        C.append([d[k][0] for k in range(deg+1)])
    # C[c][r] = coeff of binom(t2,c)*binom(t1,r)
    return C
def deg_in_each(C):
    d1=max([r for c in range(len(C)) for r in range(len(C[c])) if C[c][r]!=0], default=0)
    d2=max([c for c in range(len(C)) for r in range(len(C[c])) if C[c][r]!=0], default=0)
    return d1,d2
def peval2(C,t1,t2,m):
    tot=0; bc=1
    for c in range(len(C)):
        if c>0: bc=bc*(t2-c+1)//c
        br=1; s=0
        for r in range(len(C[c])):
            if r>0: br=br*(t1-r+1)//r
            s=(s+C[c][r]*br)%m
        tot=(tot+s*bc)%m
    return tot%m
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
    # per-wire joint pass to reach the 2-condition residue
    for outer in range(6):
        left=relift(vv); r=E.run(vv)
        stuck=[a for a in left if r[E.residx[a]]%p==0]
        if len(stuck)<=2: break
        cand=collections.defaultdict(list)
        for a in stuck:
            ws=[q for q in vars_of(E.atoms[a]) if q in SHIFT]
            ws+=[q for q in atomvalvars[a] if q in SHIFT and q not in ws]
            for w in ws:
                if influences(vv,a,w): cand[w].append(a)
        done=set()
        for w,ats in sorted(cand.items(),key=lambda kv:-len(kv[1])):
            ats=[a for a in ats if a not in done]
            if not ats: continue
            t=solve_group(vv,ats,w)
            if t is not None: vv[w]+=p*t; done.update(ats)
    left=relift(vv); r=E.run(vv)
    stuck=[a for a in left if r[E.residx[a]]%p==0]
    print('residue: %d conditions'%len(stuck),flush=True)
    # ---- component sizes of the shares-a-condition graph ----
    cand=collections.defaultdict(list)
    for a in stuck:
        ws=[q for q in vars_of(E.atoms[a]) if q in SHIFT]
        ws+=[q for q in atomvalvars[a] if q in SHIFT and q not in ws]
        for w in ws:
            if influences(vv,a,w): cand[w].append(a)
    adj=collections.defaultdict(set)
    for w,ats in cand.items():
        for x in ats:
            for y in ats: adj[x].add(y)
    seen=set(); comps=[]
    for a in stuck:
        if a in seen: continue
        st=[a]; c=set()
        while st:
            x=st.pop()
            if x in seen: continue
            seen.add(x); c.add(x); st+=list(adj[x])
        comps.append(c)
    print('COMPONENT SIZES at |S|=17:',sorted(len(c) for c in comps),flush=True)
    if len(stuck)!=2:
        print('not the 2-condition case; stopping'); sys.exit()
    a1,a2=stuck; c1=abs(SL[a1])//p; c2=abs(SL[a2])//p
    W=sorted({w for w,ats in cand.items() if ats})
    w1,w2=W[0],W[1] if len(W)>1 else W[0]
    print('atoms c1=%d %s ; c2=%d %s'%(c1,factor(c1),c2,factor(c2)),flush=True)
    print('coprime moduli? %s'%(gcd(c1,c2)==1),flush=True)
    print('wires available:',W,flush=True)
    for (x,y) in itertools.combinations(W,2):
        C1=fit2(vv,E.residx[a1],x,y); C2=fit2(vv,E.residx[a2],x,y)
        if C1 is None or C2 is None: continue
        print('  wires (x%d,x%d): deg(atom1)=%s  deg(atom2)=%s'%(x,y,deg_in_each(C1),deg_in_each(C2)),flush=True)
        # prime-by-prime: each prime constrains ONE atom (moduli coprime) -> underdetermined
        pairs1=[]; pairs2=[]; ok=True
        for (C,c,which) in ((C1,c1,1),(C2,c2,2)):
            for q,e in factor(c).items():
                m=q**e; found=None
                for t1 in range(min(m,200000)):
                    for t2 in range(m):
                        if peval2(C,t1,t2,m)==0: found=(t1,t2); break
                    if found: break
                if found is None: ok=False; break
                (pairs1 if which else pairs1).append((found[0],m))
                pairs2.append((found[1],m))
                pairs1[-1]=(found[0],m)
            if not ok: break
        if not ok: print('    no solution at some prime power'); continue
        t1=crt_list(pairs1); t2=crt_list(pairs2)
        if t1 is None or t2 is None: print('    CRT inconsistent'); continue
        v1=probe(vv,E.residx[a1],[x,y],[t1,t2]); v2=probe(vv,E.residx[a2],[x,y],[t1,t2])
        print('    t1=%d t2=%d -> atom1 %% c1*p ==0: %s   atom2 %% c2*p ==0: %s'%(
            t1,t2,v1%(c1*p)==0,v2%(c2*p)==0),flush=True)
        if v1%(c1*p)==0 and v2%(c2*p)==0:
            vv[x]+=p*t1; vv[y]+=p*t2
            l2=relift(vv); r2=E.run(vv)
            s2=[a for a in l2 if r2[E.residx[a]]%p==0]
            print('    *** BIVARIATE SOLVE SUCCEEDED: %d undischarged, %d nonzero atoms'%(
                len(s2),sum(1 for q in r2 if q)),flush=True)
            break
