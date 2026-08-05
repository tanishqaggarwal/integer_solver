"""Which free inputs move a given variable, and what else do they disturb?"""
import pickle, sys, time
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
NV=38748
checks=[a for a in range(len(polys)) if a not in atom_out]
freeinp=[x for x in range(NV) if x not in definer]
J=pickle.load(open('jac.pkl','rb'))['J']

if __name__=='__main__':
    targets=[int(x) for x in sys.argv[1:]] or [12186]
    v0=H.load_assignment('S0.json')
    res={t:[] for t in targets}
    for f in freeinp:
        v=list(v0); ch,_=ripple(v,{f:v0[f]+1})
        for t in targets:
            if t in ch: res[t].append((f, v[t]-v0[t], sorted(J.get(f,{}))))
    for t in targets:
        print(f'### free inputs moving x_{t}: {len(res[t])}')
        for f,d,ck in res[t][:40]:
            print(f'   x_{f}: dx={d}  disturbs checks {ck}')
