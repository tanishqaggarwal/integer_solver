"""Enumerate every variable adjacent to the defect and its true marginal cost
(which equations outside the 11 it would break)."""
import sys, pickle, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s9/pins')
from build import *

d=pickle.load(open('atoms.pkl','rb')); src=d['atom_src']; eq_terms=d['eq_terms']
atom2eq=pickle.load(open('atom2eq.pkl','rb'))
codes,_=H.load_equations()
FAILS=set(H.evaluate(codes,BASE))
freeset=set(x for x in range(NV) if x not in definer)

cand=set()
for e in FAILS:
    m,sq,tl=eq_terms[e]
    for c,a in tl: cand.add(a)
# also the sub-atoms of the square 37887's root
roots=pickle.load(open('roots.pkl','rb'))

# every variable occurring in those atoms
vars_=set()
for a in cand: vars_ |= set(u for m_ in polys[a] for u in m_)
print('variables occurring in the 11 failing equations\' atoms:', len(vars_))
print()
print(f'{"var":>8} {"free":>5} {"#atoms":>7} {"newEqs":>7}  atoms / equations outside FAILS')
rows=[]
for x in sorted(vars_):
    ats=var_atoms[x]
    E=set()
    for a in ats: E |= set(atom2eq.get(a,[]))
    new=sorted(E-FAILS)
    rows.append((len(new), x, x in freeset, ats, new))
rows.sort()
for n,x,fr,ats,new in rows:
    print(f'x_{x:<7} {"F" if fr else "-":>5} {len(ats):>7} {n:>7}  atoms={ats}  new={new[:8]}')
