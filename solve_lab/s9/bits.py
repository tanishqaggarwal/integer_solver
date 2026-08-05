"""Find boolean free inputs (b*(b-1)=0 style atoms) and test affineness of the free-input map."""
import pickle, random, sys
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
NV=38748
checks=[a for a in range(len(polys)) if a not in atom_out]
freeinp=[x for x in range(NV) if x not in definer]
freeset=set(freeinp)

def boolean_vars():
    """atoms of the form x^2-x (or 2x(1-x), x(x-1)) -> the var is boolean."""
    out=set()
    for a,P in enumerate(polys):
        vs=set(u for m in P for u in m)
        if len(vs)!=1: continue
        u=next(iter(vs))
        c2=P.get((u,u),0); c1=P.get((u,),0); c0=P.get((),0)
        if c2 and c0==0 and c1 == -c2:   # c2*(u^2-u)
            out.add(u)
    return out

if __name__=='__main__':
    B=boolean_vars()
    print(f'boolean-constrained vars: {len(B)}; free among them: {len(B & freeset)}')
    pickle.dump(sorted(B), open('boolvars.pkl','wb'))
    v0=H.load_assignment('S0.json')
    base={a:evalpoly(polys[a],v0) for a in checks}
    def eff(seeds):
        v=list(v0); ripple(v,seeds)
        return {a:evalpoly(polys[a],v)-base[a] for a in checks if evalpoly(polys[a],v)!=base[a]}
    random.seed(1)
    print('\naffineness test (joint vs sum of singles):')
    pairs=[(24548,14853),(24548,6418),(14853,6418),(24548,12553),(3484,6418)]
    for f,gv in pairs:
        A=eff({f:v0[f]+1}); Bq=eff({gv:v0[gv]+1}); Jt=eff({f:v0[f]+1, gv:v0[gv]+1})
        keys=set(A)|set(Bq)|set(Jt)
        ok=all(Jt.get(k,0)==A.get(k,0)+Bq.get(k,0) for k in keys)
        print(f'  x_{f} & x_{gv}: additive={ok}')
    print('\nboolean free vars and their current values:')
    bf=sorted(B & freeset)
    print(f'  {len(bf)} boolean free inputs')
    vals={b:v0[b] for b in bf}
    import collections
    print('  value histogram:', collections.Counter(vals.values()))
