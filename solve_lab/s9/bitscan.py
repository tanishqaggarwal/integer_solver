"""Flip each boolean free input; record core-activation and residual footprint."""
import pickle, time, sys, collections
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
P=2**256-2**32-977
NV=38748
roots=pickle.load(open('roots.pkl','rb'))
checks=[a for a in range(len(polys)) if a not in atom_out]
resid_poly={a:(roots[a] if a in roots else polys[a]) for a in checks}
boolv=set(pickle.load(open('boolvars.pkl','rb')))
freeinp=[x for x in range(NV) if x not in definer]
bfree=[b for b in freeinp if b in boolv]

if __name__=='__main__':
    v0=H.load_assignment('S0.json')
    base={a:evalpoly(Pp,v0) for a,Pp in resid_poly.items()}
    res={}; t0=time.time()
    for i,b in enumerate(bfree):
        v=list(v0); ch,_=ripple(v,{b: 1-v0[b]})
        touched=set()
        for u in ch: touched.update(var_atoms[u])
        nz=[]
        for a in touched:
            if a in resid_poly and evalpoly(resid_poly[a],v)!=0: nz.append(a)
        for a,x in base.items():
            if x and a not in touched: nz.append(a)
        res[b]=(v[15298], len(ch), sorted(set(nz)))
        if i%200==0: print(f'{i}/{len(bfree)} {time.time()-t0:.0f}s',file=sys.stderr)
    pickle.dump(res, open('bitscan.pkl','wb'))
    kill=[(b,r) for b,r in res.items() if r[0]==0]
    print(f'bits that switch off the core (x_15298=0): {len(kill)}')
    for b,r in sorted(kill, key=lambda t:len(t[1][2]))[:20]:
        print(f'   x_{b}: changed={r[1]} vars, {len(r[2])} nonzero residuals {r[2][:14]}')
    print()
    best=sorted(res.items(), key=lambda t: len(t[1][2]))[:15]
    print('flips with fewest residuals overall:')
    for b,r in best: print(f'   x_{b}: x_15298={r[0]} changed={r[1]} nz={len(r[2])} {r[2][:12]}')
