"""Unique-Q-solution denominators for a given enlargement."""
import sys, json, collections; sys.path.insert(0,'.')
from fractions import Fraction as F
import env, lib as L
from agrow import model
P=env.P

def qsolve(rows,nk):
    aff=[(e,c,lin) for e,c,lin,hq in rows if not hq]
    mat=[[F(lin.get(j,0)) for j in range(nk)]+[F(-c)] for e,c,lin in aff]
    nr=len(mat); r=0; piv=[]
    for col in range(nk):
        pr=None
        for i in range(r,nr):
            if mat[i][col]!=0: pr=i;break
        if pr is None: continue
        mat[r],mat[pr]=mat[pr],mat[r]
        pv=mat[r][col]; mat[r]=[x/pv for x in mat[r]]
        for i in range(nr):
            if i!=r and mat[i][col]!=0:
                f=mat[i][col]; mat[i]=[a-f*b for a,b in zip(mat[i],mat[r])]
        piv.append(col); r+=1
    incons=[aff[i][0] for i in range(r,nr) if mat[i][nk]!=0]
    sol=[None]*nk
    for i,c in enumerate(piv): sol[c]=mat[i][nk]
    free=[j for j in range(nk) if sol[j] is None]
    return sol,free,incons,r

def report(name,extra):
    A,K,R,rows,QUAD=model(extra)
    sol,free,incons,r=qsolve(rows,len(K))
    print('%-26s knobs=%d rank=%d free=%d incons=%s'%(name,len(K),r,len(free),incons[:5]))
    if incons: return
    dens=collections.Counter()
    bad=[]
    for j,u in enumerate(K):
        if sol[j] is None: continue
        d=sol[j].denominator
        if d==1: dens['1']+=1
        elif d==P: dens['p']+=1; bad.append((u,'p'))
        elif d%P==0: dens['p*%d'%(d//P)]+=1; bad.append((u,'p*%d'%(d//P)))
        else: dens[str(d)]+=1; bad.append((u,str(d)))
    print('    denominators:',dict(dens))
    if bad: print('    non-integral knobs:',bad)
    else: print('    *** ALL INTEGRAL ***')
    return sol,K

if __name__=='__main__':
    for name,ex in [('base+37887+41906',[37887,41906]),
                    ('+29426',[37887,41906,29426]),
                    ('+41972',[37887,41906,41972]),
                    ('+29426+41972',[37887,41906,29426,41972]),
                    ('+29090',[37887,41906,29090]),
                    ('+36085',[37887,41906,36085]),
                    ('+29090+36085',[37887,41906,29090,36085]),
                    ('+all',[37887,41906,29426,41972,29090,36085])]:
        report(name,ex)
