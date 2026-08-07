"""bn_base: baseline residual of the 39026 partial + structural role of boolean atoms."""
import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad

BEST = os.path.join(LAB,'best','new_instance_partial_39026.json')
v = L.load(BEST)
av = L.all_atom_values(v)
fail = L.failing_eqs(av)
nz = [a for a in range(L.NA) if av[a]]
print('base score', L.NEQ-len(fail), 'failing', fail)
print('nonzero atoms', len(nz), nz)
for a in nz:
    print(f'  a{a} val={av[a]} eqs={sorted(L.atom2eq.get(a,{}).items())[:8]} nvars={len(L.avars[a])} src={L.atom_src[a][:90]}')
print()
for e in fail:
    m,sq,co = L.eq_atoms[e]
    s = sum(c*av[a] for a,c in co.items())
    print(f'  eq{e} mult={m} sq={sq} natoms={len(co)} S={s} val={m*(s*s if sq else s)}')
    print(f'     nz members: {[(a,co[a],av[a]) for a in co if av[a]]}')

# structural: which boolean atoms are gates?
cen = json.load(open(os.path.join(HERE,'bn_census.json')))
bools = {int(a):tuple(t) for a,t in cen['bools'].items()}
gate = [a for a in bools if a in L.atom_out]
print()
print('boolean atoms that are GATES (in atom_out):', len(gate))
FREESET=set(ad.FREE)
# For each boolean atom, how many equations contain it?
h = collections.Counter(len(L.atom2eq.get(a,{})) for a in bools)
print('eqs-per-boolean-atom histogram:', sorted(h.items()))
# boolean atoms on FREE vars
bf = [a for a,(u,c) in bools.items() if u in FREESET]
print('boolean atoms on free vars:', len(bf))
h = collections.Counter(len(L.atom2eq.get(a,{})) for a in bf)
print('  eqs-per-atom histogram (free):', sorted(h.items()))
