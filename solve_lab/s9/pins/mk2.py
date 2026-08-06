"""Realise a knob-set solution directly (no ripple: the model covers every affected atom)."""
import sys, pickle, itertools, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s9/pins')
from build import *
from opt import maxzero, model, CAND
codes,_=H.load_equations()

def realise(V, sol):
    v=list(BASE)
    for var,dv in zip(sorted(V), sol): v[var]=BASE[var]+dv
    return v

def run(V, save=None):
    r=maxzero(V)
    if r is None: return None
    S,Vl,best=r
    k,sub,sol=best
    v=realise(Vl,sol)
    f=H.evaluate(codes,v)
    print(f'V={Vl}\n  |S|={len(S)} predicted zeroed={k} predicted failing={len(S)-k}')
    print(f'  ACTUAL: satisfied {len(codes)-len(f)}/{len(codes)}  ({len(f)} failing) {f}')
    if save and len(f)<11:
        H.save_assignment(v,save); print('  saved ->',save)
    return len(f), v

if __name__=='__main__':
    import ast
    V=[9413,28730]+ast.literal_eval(sys.argv[1]) if len(sys.argv)>1 else [9413,28730,642,29854,31864]
    run(V, 'pins/cand.json')
