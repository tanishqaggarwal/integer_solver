"""Construction: activate x_8599 (frees x_12186 mod p via x_5096), pin x_12186=x_14853=K1,
let the canonical gates close C1/C2, then repair the remainder."""
import pickle, sys, time
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
P=2**256-2**32-977
roots=pickle.load(open('roots.pkl','rb'))
checks=[a for a in range(len(polys)) if a not in atom_out]
resid_poly={a:(roots[a] if a in roots else polys[a]) for a in checks}
K1=33472904810391811973223207617762334363023286939839396241234196646906030803538671321618319

def nz_resid(v):
    return [a for a,Pp in resid_poly.items() if evalpoly(Pp,v)!=0]

def attempt(b, rounds=12, verbose=False):
    v=H.load_assignment('../best/new_instance_partial_39022.json')
    # 1) activate the x_8599 branch and drive x_12186 -> K1 ; 2) align x_14853 ; 3) close C1/C2
    ripple(v, {b: 1, 5096: K1, 33612: 0})
    ripple(v, {14853: v[12186]})
    ripple(v, {7068: v[2099]+7376877*v[642], 4432: v[19964]+v[28730]})
    ripple(v, {24548: v[25442]})
    ok, hist = repair_loop(v, rounds=rounds, verbose=verbose)
    return v, nz_resid(v), hist

if __name__=='__main__':
    bits=pickle.load(open('hits8599.pkl','rb'))
    if len(sys.argv)>1: bits=[int(sys.argv[1])]
    best=None; t0=time.time()
    for i,b in enumerate(bits):
        try:
            v,nz,hist=attempt(b)
        except Exception as e:
            print(f'x_{b}: ERROR {e}'); continue
        if best is None or len(nz)<len(best[1]): best=(b,nz,v)
        if len(nz)<=6 or i%20==0:
            print(f'x_{b}: {len(nz)} nonzero residuals {sorted(nz)[:10]}  hist={hist}  [{time.time()-t0:.0f}s]')
    b,nz,v=best
    print(f'\nBEST bit x_{b}: {len(nz)} residuals {sorted(nz)}')
    H.save_assignment(v, 'construct_best.json')
