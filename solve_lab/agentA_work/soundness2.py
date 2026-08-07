"""Part 2 of the soundness check: (a) confirm the RAW atom-level relaxation IS vacuous in
my parse too; (b) show explicitly that the vacuity witnesses are NOT in the image of my
knob->atom map, which is what makes my code different."""
import sys, collections; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from regsolve2 import build
P=env.P; Q=(1<<61)-1
# (a) raw relaxation
one=[a for a in range(L.NA) if len(L.atom2eq.get(a,{}))==1]
print('(a) RAW relaxation over arbitrary atom vectors:')
print('    %d atoms occur in exactly one equation, so the atom vector e_a has ||M e_a||_0 = 1.'%len(one))
print('    => min over nonzero integer atom vectors of ||M a||_0 = 1.  VACUOUS, confirmed')
print('    independently in my parse.  Example witnesses: %s'%one[:6])
priv=[a for a in one if any(len(L.var_atoms[u])==1 for u in L.avars[a])]
print('    of those %d single-equation atoms, how many carry a private variable? %d'%(len(one),len(priv)))
print()
v=L.load('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
fe=L.failing_eqs(L.all_atom_values(v))
print('(b) are the vacuity witnesses reachable in MY construction?')
for LEV in [int(x) for x in sys.argv[1:]] or [6,16]:
    A=set(a for e in fe for a in L.eq_atoms[e][2])
    for _ in range(LEV):
        R=set()
        for a in A: R|=set(L.atom2eq[a])
        A=set(a for e in R for a in L.eq_atoms[e][2])
    K,Rr,rows=build(v,A); nk=len(K)
    inwin=[a for a in one if a in A]
    # linear part of each atom in the knobs
    ki={u:i for i,u in enumerate(K)}
    def atom_lin(a):
        lin=collections.defaultdict(int)
        for m,c in L.polys[a].items():
            ks=[u for u in m if u in ki]
            if len(ks)==1:
                t=c
                for u in m:
                    if u!=ks[0]: t*=v[u]
                lin[ki[ks[0]]]+=t
        return {i:x for i,x in lin.items() if x}
    Alist=sorted(A)
    movable=[a for a in Alist if atom_lin(a)]
    print('  L=%-3d atoms=%-5d knobs=%-4d : single-equation atoms inside the window = %d'%(
        LEV,len(A),nk,len(inwin)))
    print('        atoms the knobs can move at all: %d of %d (the other %d are FROZEN at 0)'%(
        len(movable),len(A),len(A)-len(movable)))
    print('        image of the knob->atom map is a rank-<=%d sublattice of Z^%d'%(nk,len(A)))
    if inwin:
        bad=[a for a in inwin if atom_lin(a)]
        print('        single-equation atoms the knobs can move: %d  %s'%(len(bad),bad[:5]))
    else:
        print('        => NO single-equation atom is even present, so no weight-1 support')
        print('           can arise from that mechanism in this window.')
