"""Close activated load pins  bit*(x_B - HUGE) = s*x_C  by setting the free x_B := HUGE. Iterate."""
import pickle, sys, ast
import harness as H
import poly as PY
exec(open('repair.py').read().split('if __name__')[0])
P=2**256-2**32-977
roots=pickle.load(open('roots.pkl','rb'))
d=pickle.load(open('atoms.pkl','rb')); src=d['atom_src']
checks=[a for a in range(len(polys)) if a not in atom_out]
rp={a:(roots[a] if a in roots else polys[a]) for a in checks}
NV=38748
freeset=set(x for x in range(NV) if x not in definer)
K1=33472904810391811973223207617762334363023286939839396241234196646906030803538671321618319
def nz(v): return sorted(a for a,Pp in rp.items() if evalpoly(Pp,v)!=0)

def pin_fix(a, v):
    """If atom a is a pin with a free pinned var solvable exactly, return {var: value}."""
    Pp = rp[a]
    cands=[]
    for t in set(u for m in Pp for u in m):
        if t not in freeset: continue
        c=0; bad=False
        for m,cc in Pp.items():
            if len(m)==1 and m[0]==t: c+=cc
            elif t in m: bad=True
        if bad or c==0: continue
        old=v[t]; v[t]=0; rest=evalpoly(Pp,v); v[t]=old
        if rest % c: continue
        cands.append((len(var_atoms[t]), t, -rest//c))
    if not cands: return None
    cands.sort(); return {cands[0][1]: cands[0][2]}

if __name__=='__main__':
    def build(bits):
        v=H.load_assignment('../best/new_instance_partial_39022.json')
        s={b:1 for b in bits}; s[5096]=K1; s[33612]=0
        ripple(v,s); ripple(v,{14853:v[12186]})
        ripple(v,{7068:v[2099]+7376877*v[642],4432:v[19964]+v[28730]}); ripple(v,{24548:v[25442]})
        return v
    v=build([2527,1502])
    hist=[]
    for it in range(25):
        cur=nz(v); hist.append(len(cur))
        print(f'iter {it}: {len(cur)} residuals {cur[:12]}')
        if not cur: print('ALL CLEAR'); break
        seeds={}
        for a in cur:
            f=pin_fix(a,v)
            if f:
                for k,val in f.items():
                    if k not in seeds: seeds[k]=val
        if not seeds: print('no pin-style repair available -> stuck'); break
        ripple(v, seeds)
        ripple(v,{7068:v[2099]+7376877*v[642],4432:v[19964]+v[28730]})
    print('history:', hist)
    codes,_=H.load_equations(); f=H.evaluate(codes,v)
    print(f'EQUATIONS: {len(codes)-len(f)}/{len(codes)} ({len(f)} failing)')
