"""Two-bit construction: activate x_8599 branch AND zero the second core's controls."""
import pickle, sys, time, itertools
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
P=2**256-2**32-977
roots=pickle.load(open('roots.pkl','rb'))
checks=[a for a in range(len(polys)) if a not in atom_out]
rp={a:(roots[a] if a in roots else polys[a]) for a in checks}
K1=33472904810391811973223207617762334363023286939839396241234196646906030803538671321618319
def nz(v): return sorted(a for a,Pp in rp.items() if evalpoly(Pp,v)!=0)

def build(bits, k1=K1):
    v=H.load_assignment('../best/new_instance_partial_39022.json')
    seeds={b:1 for b in bits}; seeds[5096]=k1; seeds[33612]=0
    ripple(v, seeds)
    ripple(v, {14853: v[12186]})
    ripple(v, {7068: v[2099]+7376877*v[642], 4432: v[19964]+v[28730]})
    ripple(v, {24548: v[25442]})
    return v

if __name__=='__main__':
    A=pickle.load(open('hits8599.pkl','rb'))
    Bcand=[47,112,1438,1502,1544,1571,1746,1931,2055,3106,3216,3414,3707,4381]
    print('intersection hits8599 & B:', sorted(set(A)&set(Bcand)))
    best=None; t0=time.time(); results=[]
    for b1 in A:
        for b2 in Bcand:
            if b1==b2: continue
            v=build([b1,b2])
            n=nz(v)
            results.append((len(n), b1, b2, n))
            if best is None or len(n)<best[0]: best=(len(n),b1,b2,n,v)
    results.sort()
    print(f'{len(results)} pairs tried in {time.time()-t0:.0f}s')
    for r in results[:12]: print(f'   ({r[1]},{r[2]}): {r[0]} residuals {r[3][:10]}')
    n,b1,b2,nn,v = best
    print(f'\nBEST pair (x_{b1}, x_{b2}): {n} residuals {nn}')
    H.save_assignment(v,'construct3_best.json')
    codes,_=H.load_equations(); f=H.evaluate(codes,v)
    print(f'EQUATIONS: {len(codes)-len(f)}/{len(codes)} ({len(f)} failing)')
