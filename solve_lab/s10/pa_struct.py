import os, sys, collections, json
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L
# equation size histogram
h=collections.Counter(len(co) for m,sq,co in L.eq_atoms)
print('eq size hist:', sorted(h.items())[:20])
print('n eqs with <=3 atoms:', sum(c for k,c in h.items() if k<=3))
print('n eqs with 1 atom:', h.get(1,0), ' 2 atoms:', h.get(2,0), ' 3:', h.get(3,0))
# squares
print('is_square count:', sum(1 for m,sq,co in L.eq_atoms if sq))
# footprint-1 atoms: which equations do they land in?
f1=[a for a in range(L.NA) if len(L.atom2eq[a])==1]
print('footprint-1 atoms:', len(f1))
eqc=collections.Counter(next(iter(L.atom2eq[a])) for a in f1)
print('eqs hosting footprint1 atoms:', len(eqc), 'hist of #f1 per eq:', sorted(collections.Counter(eqc.values()).items()))
# print a few
for a in f1[:6]:
    e=next(iter(L.atom2eq[a])); m,sq,co=L.eq_atoms[e]
    print(f'  a{a} (gate={a in L.atom_out}) eq{e} size={len(co)} src={L.atom_src[a][:90]}')
print()
for r in [(35750,35770),(5760,5780),(22225,22240),(19085,19095),(10930,10942)]:
    for a in range(*r):
        if a<L.NA:
            print(f'a{a} neq={len(L.atom2eq[a]):<3} gate={1 if a in L.atom_out else 0} {L.atom_src[a][:110]}')
    print()
