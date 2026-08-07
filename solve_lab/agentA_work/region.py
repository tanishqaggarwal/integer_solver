"""Exact region model around the seven residual atoms.
Level-k closure: atoms A, equations R, knobs = vars whose every atom lies in A."""
import sys, collections, json; sys.path.insert(0,'.')
import env, lib as L
P=env.P
v0=env.load_best(); av0=L.all_atom_values(v0)

def close(A):
    A=set(A)
    knobs=[u for u in set(u for a in A for u in L.avars[a]) if all(x in A for x in L.var_atoms[u])]
    R=sorted(set(e for a in A for e in L.atom2eq[a]))
    return sorted(A), sorted(knobs), R

def report(A,label):
    A,knobs,R=close(A)
    # which atoms appear in R but are outside A?
    outside=set()
    for e in R: outside |= set(L.eq_atoms[e][2])-set(A)
    nzout=[a for a in outside if av0[a]]
    print('%-22s atoms=%-5d knobs=%-4d eqs=%-5d foreign_atoms_in_eqs=%-5d (nonzero %d)'%(
        label,len(A),len(knobs),len(R),len(outside),len(nzout)))
    return A,knobs,R,outside

E=sorted(set(e for a in env.SEVEN for e in L.atom2eq[a]))
A0=sorted(set(a for e in E for a in L.eq_atoms[e][2]))
lvl=A0
report(lvl,'L0 (region of E)')
for k in range(1,7):
    # grow: add every atom touched by a variable of the current atom set
    new=set(lvl)
    for a in lvl: new |= set(x for u in L.avars[a] for x in L.var_atoms[u])
    lvl=sorted(new)
    A,knobs,R,outside=report(lvl,'L%d (var-closure)'%k)
    if len(R)>4000: break
