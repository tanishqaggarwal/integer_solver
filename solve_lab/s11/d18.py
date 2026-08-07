import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
v=load_raw(os.path.join(HERE,'data','finish3_named.json'))
FREE=set(u for u in range(L.NVARS) if u not in L.definer)
def sh(x):
    s=str(x); return s if len(s)<18 else s[:7]+'..'+s[-4:]+f'<{len(s)}d>'
BR=[26719,26721,26723]
EQS=sorted(set().union(*[set(L.atom2eq.get(a,{})) for a in BR]))
ATOMS=sorted(set().union(*[set(L.eq_atoms[e][2]) for e in EQS]))
print("broken:",BR," equations touched:",len(EQS)," atoms in them:",len(ATOMS))
for a in ATOMS:
    Pp=L.polys[a]
    t=' + '.join(f"{c}*{'*'.join('x%d'%z for z in m)}" for m,c in Pp.items())
    print(f"a{a} out={str(L.atom_out.get(a)):12s} val={sh(L.evalpoly(Pp,v))} neq={len(L.atom2eq.get(a,{}))}  {t[:150]}")
print()
VARS=sorted(set().union(*[set(L.avars[a]) for a in ATOMS]))
print(f"{len(VARS)} variables; fully local ones (no atom outside this set):")
for u in VARS:
    out=[a for a in L.var_atoms[u] if a not in set(ATOMS)]
    if not out:
        print(f"   x{u} {'FREE' if u in FREE else 'DEF a%d'%L.definer[u]} val={sh(v[u])} atoms={L.var_atoms[u]}")
