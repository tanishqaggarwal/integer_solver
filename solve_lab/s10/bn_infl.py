"""bn_infl: can any boolean variable influence the 7 residual atoms / 7 failing eqs?

Computes the ancestor closure (in the gate DAG) of the variables appearing in the
7 nonzero atoms of the 39026 state, and intersects with boolean variables.
Also audits the two candidate blocks quoted in the task brief.
"""
import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'
sys.path.insert(0, HERE)
import bn_lib as B
import lib as L, tools as T, ad

bools=B.bools_map()
bvar={u for u,_ in bools.values()}
NZ=[22229,22230,35758,35759,35760,35761,35762]
FAIL=[12231,12270,12350,14584,18673,22044,29125]

seed=set()
for a in NZ: seed |= set(L.avars[a])
print('vars in the 7 residual atoms:',sorted(seed))

# ancestor closure: x_t <- vars of definer[t]
anc=set(seed); st=list(seed)
while st:
    u=st.pop()
    d=L.definer.get(u)
    if d is None: continue
    for w in L.avars[d]:
        if w not in anc: anc.add(w); st.append(w)
print('ancestor closure size:',len(anc))
print('boolean vars among ancestors:',len(anc & bvar))
print('FREE boolean vars among ancestors:',len(anc & bvar & B.FREESET))
print('FREE vars among ancestors:',len(anc & B.FREESET))

# also: atoms of the failing equations, and their ancestor closure
seed2=set()
for e in FAIL:
    m,sq,co=L.eq_atoms[e]
    for a in co: seed2 |= set(L.avars[a])
anc2=set(seed2); st=list(seed2)
while st:
    u=st.pop()
    d=L.definer.get(u)
    if d is None: continue
    for w in L.avars[d]:
        if w not in anc2: anc2.add(w); st.append(w)
print()
print('ancestor closure of ALL atoms in the 7 failing eqs:',len(anc2))
print('  boolean vars there:',len(anc2 & bvar),' free boolean:',len(anc2 & bvar & B.FREESET))

# audit the blocks quoted in the brief
print()
print('=== audit of quoted candidate blocks ===')
for a in (33516,33517,27821,27822,33522):
    isb = a in bools
    print(f'  a{a}: boolean={isb} vars={sorted(L.avars[a])} '
          f'|E|={len(L.atom2eq.get(a,{}))} src={L.atom_src[a][:70]}')
for x in (29570,27026,33095,24267):
    print(f'  x_{x}: boolean_var={x in bvar} free={x in B.FREESET} natoms={len(L.var_atoms[x])}')

json.dump({'anc_resid':sorted(anc),'bool_anc':sorted(anc&bvar),
           'freebool_anc':sorted(anc&bvar&B.FREESET),
           'bool_anc2':sorted(anc2&bvar&B.FREESET)},
          open(os.path.join(HERE,'bn_infl.json'),'w'))
