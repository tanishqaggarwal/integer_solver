"""bl_cancel: satisfy a failing equation by CANCELLATION rather than by zeroing
the residual atoms.

Each failing equation is a small-integer combination of atoms.  Most of the atoms
in it are currently zero gate atoms.  If some atom a appears ONLY in equations
that already fail, we may give it any value we like; choosing it to cancel one
equation gains +1 with no loss -- provided the variable that realises it is not
used anywhere else.
"""
import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from bl_core import Frame, CANON, F2, pot
P = 2**256-2**32-977

w = L.load(os.path.join(LAB,'best','new_instance_partial_39026.json')); F2.fwd(w)
av = L.all_atom_values(w)
FAIL = set(L.failing_eqs(av))
print(f'failing equations: {sorted(FAIL)}')
EQA = sorted(set().union(*[set(L.eq_atoms[i][2]) for i in FAIL]))
print(f'atoms appearing in them: {len(EQA)}  {EQA}')
print(f'\n{"atom":>7} {"value":>6} {"#eqs":>5} {"#failing":>8} {"#satisfied":>10}  out       vars')
free_atoms = []
for a in EQA:
    eqs = set(L.atom2eq.get(a, ()))
    nf = len(eqs & FAIL); ns = len(eqs - FAIL)
    print(f'  a{a:<6} {"NZ" if av[a] else "0":>4} {len(eqs):>5} {nf:>8} {ns:>10}  '
          f'{str(L.atom_out.get(a)):<10}{sorted(L.avars[a])}')
    if ns == 0 and not av[a]: free_atoms.append(a)
print(f'\natoms that occur ONLY in already-failing equations: {free_atoms}')

# for each such atom, and each failing equation, what value would satisfy it?
for a in free_atoms:
    t = L.definer.get(a) if False else None
    o = L.atom_out.get(a)
    print(f'\n--- a{a}: out={o} src={L.atom_src[a][:90]}')
    if o:
        tv = o[1]
        others = [b for b in L.var_atoms[tv] if b != a]
        oeqs = set()
        for b in others: oeqs |= set(L.atom2eq.get(b, ()))
        print(f'    its output var x_{tv} also appears in atoms {others} '
              f'-> {len(oeqs)} equations, {len(oeqs - FAIL)} of them currently satisfied')
    for i in sorted(FAIL):
        m, sq, co = L.eq_atoms[i]
        if a not in co: continue
        s = sum(c * av[b] for b, c in co.items())
        c = co[a]
        need = -s + c * av[a]
        if need % c == 0:
            print(f'    eq{i}: a{a} := {str(need // c)[:44]}... would satisfy it '
                  f'(coef {c}, exact)')
        else:
            print(f'    eq{i}: coef {c} does not divide the residue -- not solvable alone')
