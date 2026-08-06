"""Anatomy of the 11 failing equations: which atoms compose them, and which atoms are
'private' (appear in few equations) so could absorb the residual."""
import sys, pickle, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s9/pins')
from build import *

d=pickle.load(open('atoms.pkl','rb')); src=d['atom_src']; eq_terms=d['eq_terms']
atom2eq=pickle.load(open('atom2eq.pkl','rb'))
codes,_=H.load_equations()
FAILS=H.evaluate(codes,BASE)
print('failing equations:',FAILS)
roots=pickle.load(open('roots.pkl','rb'))
freeset=set(x for x in range(NV) if x not in definer)

for e in FAILS:
    terms=eq_terms[e]
    r=H.resid(codes,BASE,e)
    print(f'\n=== eq {e}: residual={r}  ({len(terms)} atom terms)')
    priv=[]
    for coef,a in terms:
        mult=len(atom2eq.get(a,[]))
        val=evalpoly(polys[a],BASE)
        gate = a in atom_out
        if mult<=2:
            priv.append((a,coef,mult,val,gate))
    print(f'   atoms with equation-multiplicity <=2: {len(priv)}')
    for a,coef,mult,val,gate in priv[:12]:
        issq = a in roots
        print(f'      atom {a} coef={coef} mult={mult} val0={val} gate={gate} square={issq}  {src[a][:90]}')
