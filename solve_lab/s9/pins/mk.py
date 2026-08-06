"""Realise a (z,t,h) solution as an assignment and score it."""
import sys, pickle, itertools, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s9/pins')
from build import *

codes,_=H.load_equations()

def realise(z, t, h):
    v=list(BASE)
    v[28730]=z
    v[9413]=h
    if t:
        v[7068]=BASE[7068]+P*t
        ripple(v,{})           # no-op; do the mirror by hand below
        # x_1308 tracks x_7068 ; let the canonical ripple propagate from x_7068
        ch,_=ripple_from(v,7068)
    return v

def ripple_from(v, u):
    return ripple(v, {u: v[u]})

if __name__=='__main__':
    res=pickle.load(open('pins/dioph2.pkl','rb'))
    best=None
    for sub,s in res[1]:
        z,t,h=s
        v=list(BASE); v[28730]=z; v[9413]=h
        if t:
            v[7068]=BASE[7068]+P*t
            ripple(v,{7068:BASE[7068]+P*t})
        f=H.evaluate(codes,v)
        A=evalpoly(polys[22229],v); C=evalpoly(polys[22230],v); B=evalpoly(polys[22231],v)
        print(f'target {sub}: fails={len(f)} {f}   A0={A==0} C0={C==0} B0={B==0}')
        if best is None or len(f)<len(best[0]): best=(f,v,sub)
    f,v,sub=best
    print(f'\nBEST: {len(codes)-len(f)}/{len(codes)}  ({len(f)} failing) via {sub}')
    H.save_assignment(v,'pins/best.json')
