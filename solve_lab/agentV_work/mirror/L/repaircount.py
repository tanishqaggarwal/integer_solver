"""How many of the 927 c>1 integer conditions are actually exercised and discharged?"""
import sys, pickle, json, collections
from math import gcd
src=open('/home/user/integer_solver/solve_lab/agentV_work/mirror/L/mkassign2.py').read().split('#MAINSTART')[0]
exec(src)
CGT={a for a,s in SL.items() if s and s%p==0 and abs(s)//p>1}
print('atoms carrying a c>1 integer condition: %d'%len(CGT))
import random
rnd=random.Random(7)
for S in ([M['live'][0]], [24601,2081], rnd.sample(M['live'],17), rnd.sample(M['live'],128)):
    v,isl,valn=assignment(set(S),ORIENT); v[24468]=T1; v[18956]=T2
    vv=[0]*NV
    for k,x in v.items(): vv[k]=x
    tot=0; rounds=0; hit=set()
    for rd in range(60):
        bad=relift(vv)
        rounds=rd+1
        realbad=[a for a in bad if E.run(vv)[E.residx[a]]%p==0]
        if not realbad: break
        r=E.run(vv); fixed=0
        for a in realbad:
            i=E.residx[a]; cur=r[i]; sm=abs(SL[a])
            hit.add(a)
            imm=[q for q in vars_of(E.atoms[a]) if q in SHIFT]
            for w in imm+[q for q in atomvalvars[a] if q in SHIFT and q not in imm]:
                old=vv[w]; vv[w]=old+p; d=E.run(vv)[i]-cur; vv[w]=old
                if d==0: continue
                g=gcd(d,sm)
                if cur%g: continue
                mm=sm//g
                t=(-(cur//g))*pow((d//g)%mm,-1,mm)%mm if mm>1 else 0
                vv[w]=old+p*t; fixed+=1; tot+=1; break
        if fixed==0: break
    left=relift(vv); r=E.run(vv)
    nz=[i for i,x in enumerate(r) if x]
    stillbad=[a for a in left if r[E.residx[a]]%p==0]
    print('|S|=%-4d  c>1 conditions violated & repaired: %d (over %d rounds, %d distinct atoms, all in the c>1 set: %s)'%(
        len(S),tot,rounds,len(hit),hit<=CGT))
    print('          UNDISCHARGED integer conditions remaining: %d   nonzero atoms: %d'%(len(stillbad),len(nz)))
