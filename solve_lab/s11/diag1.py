"""Global picture: gates vs checks at the checkpoint, and mod-p status."""
import sys, os, json, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw, atomval
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
LAB = os.path.join(HERE, '..')
src = sys.argv[1] if len(sys.argv)>1 else os.path.join(LAB,'best','new_instance_partial_39026.json')
v = load_raw(src)
AV = [atomval(a, v) for a in range(L.NA)]
GATE = [a for a in range(L.NA) if L.atom_out.get(a) is not None]
CHK  = [a for a in range(L.NA) if L.atom_out.get(a) is None]
print("atoms", L.NA, "gate-atoms", len(GATE), "check-atoms", len(CHK))
bg = [a for a in GATE if AV[a]!=0]
bc = [a for a in CHK  if AV[a]!=0]
print("broken gate atoms:", len(bg), bg)
print("broken check atoms:", len(bc), bc[:50])
def eqv(e):
    m,sq,co = L.eq_atoms[e]
    s = sum(c*AV[a] for a,c in co.items())
    return s
F = [e for e in range(L.NEQ) if eqv(e)!=0]
print("failing eqs", len(F), F)
for e in F:
    m,sq,co = L.eq_atoms[e]
    s = eqv(e)
    print(f"  eq{e} mult={m} sq={sq} atoms={list(co.items())} inner={s} | inner%p={s%P} | inner/p={'-' if s%P else s//P}")
# free vars
defined = set(L.definer)
free = [u for u in range(L.NVARS) if u not in defined]
print("defined vars", len(defined), "free vars", len(free))
