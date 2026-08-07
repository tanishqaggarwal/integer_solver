"""Generic residual-region exact solver.
For a state file: region = atoms in the failing equations; knobs = vars all of whose
atoms lie in the region; equations = all equations touching region atoms.
Report rank/consistency/denominators of the unique (or general) rational solution."""
import sys, json, collections; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
from fractions import Fraction as F
import env, lib as L
P=env.P

def knobpoly(a,K,v):
    ki={u:i for i,u in enumerate(K)}
    out=collections.defaultdict(int)
    for m,c in L.polys[a].items():
        t=c; mono=[]
        for u in m:
            if u in ki: mono.append(ki[u])
            else: t*=v[u]
        if t: out[tuple(sorted(mono))]+=t
    return {k:c for k,c in out.items() if c}

def analyse(v, A, verbose=True):
    av=L.all_atom_values(v)
    A=set(A); 
    K=sorted(u for u in set(u for a in A for u in L.avars[a]) if all(x in A for x in L.var_atoms[u]))
    R=sorted(set(e for a in A for e in L.atom2eq[a]))
    kp={a:knobpoly(a,K,v) for a in A}
    QUAD=set(a for a in A if any(len(m)>1 for m in kp[a]))
    def aff(a):
        Pp=kp[a]; return Pp.get((),0), {m[0]:c for m,c in Pp.items() if len(m)==1}
    rows=[]
    for e in R:
        mm,sq,co=L.eq_atoms[e]
        c=0; lin=collections.defaultdict(int); hq=False
        for a,cc in co.items():
            if a not in A:
                if av[a]!=0: hq=True
                continue
            if a in QUAD: hq=True; continue
            c0,l0=aff(a); c+=cc*c0
            for i,x in l0.items(): lin[i]+=cc*x
        rows.append((e,c,dict(lin),hq))
    return K,R,rows,QUAD

def qsolve(rows,nk):
    a=[(e,c,lin) for e,c,lin,hq in rows if not hq]
    mat=[[F(lin.get(j,0)) for j in range(nk)]+[F(-c)] for e,c,lin in a]
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
                f=mat[i][col]; mat[i]=[x-f*y for x,y in zip(mat[i],mat[r])]
        piv.append(col); r+=1
    incons=[a[i][0] for i in range(r,nr) if mat[i][nk]!=0]
    sol=[None]*nk
    for i,c in enumerate(piv): sol[c]=mat[i][nk]
    return sol,[j for j in range(nk) if sol[j] is None],incons,r,[x[0] for x in a]

def run(path,grow=0):
    v=L.load(path); av=L.all_atom_values(v); fe=L.failing_eqs(av)
    nz=[a for a in range(L.NA) if av[a]]
    print('%s  score=%d  failing=%d  nonzero atoms=%s'%(path.split('/')[-1],L.NEQ-len(fe),len(fe),nz))
    A=set(a for e in fe for a in L.eq_atoms[e][2])
    for g in range(grow+1):
        K,R,rows,QUAD=analyse(v,A)
        sol,free,incons,r,eqs=qsolve(rows,len(K))
        nq=sum(1 for x in rows if x[3])
        bad=[]
        if not incons:
            for j,u in enumerate(K):
                if sol[j] is None: continue
                d=sol[j].denominator
                if d!=1: bad.append((u,'p' if d==P else ('p*%d'%(d//P) if d%P==0 else str(d))))
        print('  g%d atoms=%-5d knobs=%-4d eqs=%-5d quadrows=%-4d rank=%-4d free=%-4d incons=%-3d nonint=%-3d %s'%(
            g,len(A),len(K),len(R),nq,r,len(free),len(incons),len(bad),bad[:6]),flush=True)
        if not incons and not bad and not free and nq==0:
            print('  *** INTEGRAL FULL SOLUTION OF THE REGION ***')
            return K,sol
        if g<grow:
            A |= set(x for a in list(A) for u in L.avars[a] for x in L.var_atoms[u])
    return None

if __name__=='__main__':
    import glob,os
    for p in sys.argv[1:]:
        try: run(p, grow=1)
        except Exception as ex: print(p,'ERR',ex)
