"""Simultaneous CRT solve over the shift parameters -- replaces the greedy round-robin.

Method inherited from P:  never brute-force over lcm(c_k) -- factor, go prime-by-prime, CRT.
Guard inherited from P:   never trust the expansion -- verify by DIRECT RECOMPUTATION.

Problem.  After the mod-p fold every atom is == 0 mod p.  For the c>1 atoms we need
c_a*p | R_a.  Knob: shift value wire w by p*t_w.  Writing R_a = p*r_a, the condition is
    r_a + (dR_a/p)(t) == 0   (mod c_a)
which is a POLYNOMIAL in t, not linear (a product of two shifted wires contributes p*t_w*t_v
after dividing by p).  So we measure the shape rather than assuming it.
"""
import sys, pickle, collections, json, itertools
from math import gcd
src=open('/home/user/integer_solver/solve_lab/agentV_work/mirror/L/mkassign2.py').read().split('#MAINSTART')[0]
exec(src)
CGT={a for a,s in SL.items() if s and s%p==0 and abs(s)//p>1}
def factor(n):
    f={}; d=2
    while d*d<=n:
        while n%d==0: f[d]=f.get(d,0)+1; n//=d
        d+=1
    if n>1: f[n]=f.get(n,0)+1
    return f
def probe(vv,a,ws,ts):
    """exact recomputation of atom a with wire shifts ws->ts (multiples of p)"""
    old=[vv[w] for w in ws]
    for w,t in zip(ws,ts): vv[w]+=p*t
    r=E.run(vv)[E.residx[a]]
    for w,o in zip(ws,old): vv[w]=o
    return r
def analyse(vv,bad):
    """for each violated atom: its c, its influencing wires, and whether it is LINEAR in them"""
    info=[]
    for a in bad:
        i=E.residx[a]; c=abs(SL[a])//p
        cur=E.run(vv)[i]
        ws=[q for q in vars_of(E.atoms[a]) if q in SHIFT]
        ws+=[q for q in atomvalvars[a] if q in SHIFT and q not in ws]
        eff=[]
        for w in ws:
            d1=probe(vv,a,[w],[1])-cur
            d2=probe(vv,a,[w],[2])-cur
            if d1==0 and d2==0: continue
            lin = (d2==2*d1)
            eff.append((w,d1,lin))
        info.append(dict(atom=a,c=c,cur=cur,r=cur//p,eff=eff,
                         all_linear=all(l for _,_,l in eff)))
    return info
if __name__=='__main__':
    import random
    rnd=random.Random(7)
    for S in ([24601,2081], rnd.sample(M['live'],17)):
        v,isl,valn=assignment(set(S),ORIENT); v[24468]=T1; v[18956]=T2
        vv=[0]*NV
        for k,x in v.items(): vv[k]=x
        # run the OLD greedy repair to reach its stuck point
        for rd in range(60):
            bad=relift(vv)
            if not bad: break
            r=E.run(vv); fixed=0
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
                    vv[w]=old+p*t; fixed+=1; break
            if fixed==0: break
        left=relift(vv); r=E.run(vv)
        stuck=[a for a in left if r[E.residx[a]]%p==0]
        print('\n|S|=%-3d greedy leaves %d undischarged'%(len(S),len(stuck)),flush=True)
        info=analyse(vv,stuck)
        for d in info:
            print('   c=%-10d r mod c = %-10d linear=%-5s  influencing wires: %s'%(
                d['c'], d['r']%d['c'], d['all_linear'],
                [(w,'d/p=%s'%(dd//p if dd%p==0 else dd)) for w,dd,_ in d['eff']][:6]),flush=True)
            print('       c factors: %s'%factor(d['c']))
