"""Exact 2-D reachable lattice for (a20215, a28647) subject to all other atoms exactly zero."""
import sys, json, collections, pickle, time, math
sys.path.insert(0,'.')
import common as C, lat2, lattice as L
import harness as H, engine as E, fast, intsolve
from flint import fmpz_mat
P=C.P
TGT=[20215,28647]

def analyse(seed,label,extra=()):
    v0,bad0,aff,atoms,hs=lat2.system(seed,extra)
    knobs=sorted(aff)
    other=[a for a in atoms if a not in TGT]
    tg=[a for a in TGT if a in atoms]
    D=lambda rows:[[aff[f].get(a,0) for f in knobs] for a in rows]
    A=D(other); b=[-bad0.get(a,0) for a in other]
    t0=time.time()
    n0,ker=intsolve.solve_int(A,b)
    print(f"[{label}] knobs={len(knobs)} other-rows={len(other)} target-rows={tg}: "
          f"{'FEASIBLE' if n0 is not None else 'INFEASIBLE'} kernel-dim={len(ker)} ({time.time()-t0:.0f}s)",flush=True)
    if n0 is None: return None
    Dt=D(tg)
    # residual on target rows at n0
    r=[bad0.get(a,0)+sum(Dt[i][j]*n0[j] for j in range(len(knobs))) for i,a in enumerate(tg)]
    # lattice generators: Dt * k for each kernel vector k
    gens=[[sum(Dt[i][j]*k[j] for j in range(len(knobs))) for i in range(len(tg))] for k in ker]
    gens=[g for g in gens if any(g)]
    print(f"   residual at particular solution: {[ (str(x)[:20]+'..%dbits'%x.bit_length()) if x else '0' for x in r]}")
    print(f"   {len(gens)} nonzero lattice generators on {tg}")
    for g in gens[:12]:
        print("      gen:",[('0' if x==0 else ('p' if abs(x)==P else ('%d*p'%(x//P) if x%P==0 else '%dbits'%x.bit_length()))) for x in g])
    if not gens:
        print("   LATTICE IS TRIVIAL -> need residual exactly 0:", all(x==0 for x in r)); return None
    M=fmpz_mat([[g[i] for g in gens] for i in range(len(tg))])   # rows=targets, cols=gens
    # membership test: solve M y = -r over Z
    Mi=[[int(M[i,j]) for j in range(M.ncols())] for i in range(M.nrows())]
    y,k2=intsolve.solve_int(Mi,[-x for x in r])
    print(f"   MEMBERSHIP of -residual in reachable lattice: {'YES -> FULL SOLVE' if y is not None else 'NO'}")
    if y is None:
        # what IS the lattice? HNF of the 2 x m generator matrix
        Hh=fmpz_mat(Mi).hnf().tolist()
        print("   HNF of generator matrix (reachable lattice basis):")
        for row in Hh[:2]:
            print("     ",[('0' if int(x)==0 else ('p' if abs(int(x))==P else ('%d*p'%(int(x)//P) if int(x)%P==0 else '%dbits/%s'%(int(x).bit_length(),'gcdp=%d'%math.gcd(abs(int(x)),P))))) for x in row[:6]])
        # index/quotient check: is -r in the lattice mod p?
        for i,a in enumerate(tg):
            print(f"   row a{a}: residual mod p = {r[i]%P}")
    return y,ker,n0,knobs,aff,v0,bad0,atoms

if __name__=='__main__':
    analyse(dict(C.BASE),'cfg0')
