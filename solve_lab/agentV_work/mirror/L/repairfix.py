"""Ordered divisibility repair: shift wires bottom-up so a fix never re-breaks a lower atom.

Each repair shifts ONE value wire w by a multiple of p.  w feeds only atoms at or ABOVE its own
position in the tree, so if atoms are processed in increasing tree height of the wire they shift,
one pass suffices and the round-robin cycling disappears.
"""
import sys, pickle, collections, json
from math import gcd
src=open('/home/user/integer_solver/solve_lab/agentV_work/mirror/L/mkassign2.py').read().split('#MAINSTART')[0]
exec(src)
NODE=M['NODE']; OUT=M['OUT']; ROOT=M['ROOT']; tree=M['tree']
parent={}
for n in NODE:
    for s,ch in (('va',NODE[n]['a']),('vb',NODE[n]['b'])): parent[ch]=n
# height of every wire = tree depth of the node that owns it (leaves deepest)
depth={ROOT:0}
def sd(n,d):
    depth[n]=d
    if tree[n] is not None:
        for c in tree[n]: sd(c,d+1)
sd(ROOT,0)
wdepth={}
for n in NODE:
    for d in OUT[n]:
        for k in ('va','vb','vab','out'): wdepth[d[k]]=depth[n]
for L in M['live']+M['dead']: wdepth.setdefault(L,depth[L])
CGT={a for a,s in SL.items() if s and s%p==0 and abs(s)//p>1}
def repair(vv,rounds=8,ordered=True):
    """returns (n_shifts, n_rounds, undischarged)"""
    tot=0; rd=0
    for rd in range(1,rounds+1):
        bad=relift(vv)
        r=E.run(vv)
        real=[a for a in bad if r[E.residx[a]]%p==0]
        if not real: break
        # order atoms by the DEPTH of the deepest shiftable wire they own: deepest first
        def key(a):
            ws=[q for q in vars_of(E.atoms[a]) if q in SHIFT]
            return -max((wdepth.get(q,0) for q in ws), default=0)
        if ordered: real.sort(key=key)
        fixed=0
        for a in real:
            i=E.residx[a]; cur=E.run(vv)[i] if ordered else r[i]
            if cur==0 or cur%p: continue
            sm=abs(SL[a])
            if cur%sm==0: continue
            imm=[q for q in vars_of(E.atoms[a]) if q in SHIFT]
            imm.sort(key=lambda q:-wdepth.get(q,0))
            for w in imm+[q for q in atomvalvars[a] if q in SHIFT and q not in imm]:
                old=vv[w]; vv[w]=old+p; d2=E.run(vv)[i]-cur; vv[w]=old
                if d2==0: continue
                g=gcd(d2,sm)
                if cur%g: continue
                mm=sm//g
                t=(-(cur//g))*pow((d2//g)%mm,-1,mm)%mm if mm>1 else 0
                vv[w]=old+p*t; fixed+=1; tot+=1; break
        if fixed==0: break
    left=relift(vv); r=E.run(vv)
    und=[a for a in left if r[E.residx[a]]%p==0]
    return tot,rd,und
if __name__=='__main__':
    import random
    rnd=random.Random(7)
    for S in ([M['live'][0]],[24601,2081],rnd.sample(M['live'],17),rnd.sample(M['live'],40)):
        for ordered in (False,True):
            v,isl,valn=assignment(set(S),ORIENT); v[24468]=T1; v[18956]=T2
            vv=[0]*NV
            for k,x in v.items(): vv[k]=x
            tot,rd,und=repair(vv,rounds=8,ordered=ordered)
            r=E.run(vv); nz=sum(1 for x in r if x)
            print('|S|=%-4d ordered=%-5s shifts=%-4d rounds=%-2d UNDISCHARGED=%-3d nonzero atoms=%d'%(
                len(S),ordered,tot,rd,len(und),nz),flush=True)
