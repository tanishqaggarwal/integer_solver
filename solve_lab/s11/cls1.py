"""Sweep all 15 (U-side, V-side) classes.

U = OR(A-leaf, B-leaf), V = OR(C-leaf, D-leaf); the mirrors are x38170 = A-leaf AND B-leaf and
x3896 = C-leaf AND D-leaf.  So the configuration space is (which of A,B fire) x (which of C,D
fire) = 16 classes, one forbidden by the OR gate.  Build a representative message for each and
score it exactly.  With the mirror OFF the checks it gates become vacuous, so this is where the
obstruction set itself changes shape.
"""
import sys, os, json, itertools, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp, forwardp
from bits5 import msg, fails, real, tree
P=L.P; sys.set_int_max_str_digits(400000)
byt=collections.defaultdict(list)
for b in real: byt[tree.get(b)].append(b)
reps={'A':byt['A'][0],'B':byt['B'][0],'C':byt['C'][0],'D':byt['D'][0]}
print("representatives:", reps)
rows=[]
for us in [(), ('A',), ('B',), ('A','B')]:
    for vs in [(), ('C',), ('D',), ('C','D')]:
        if not us and not vs: continue
        S={reps[t] for t in us+vs}
        v=msg(S); F=fails(v)
        rows.append((len(F), us, vs, v[7715], v[34554], v[38170], v[3896],
                     v[15298], v[5647], v[34606], F[:6]))
rows.sort()
print(f"{'U-side':10s} {'V-side':10s}  U V  m2 m1  chan(UV/(1-U)V/U(1-V))  failing")
for n,us,vs,U,V,m2,m1,c1,c2,c3,F in rows:
    ch = 'UV' if c1 else ('(1-U)V' if c2 else ('U(1-V)' if c3 else 'none'))
    print(f"{str(us):10s} {str(vs):10s}  {U} {V}   {m2}  {m1}   {ch:8s}  {n:3d}   {F}")
