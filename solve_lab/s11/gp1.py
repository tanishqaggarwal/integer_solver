"""GLOBAL mod-p census.

220 'wire' variables hold the value p.  Mod p they are ZERO, so every monomial containing one
dies.  Reducing the whole instance mod p therefore collapses all the quotient-witness machinery
and leaves the bare computation.  How big is what's left?
"""
import sys, os, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw, atomval
P = L.P
HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
v = load_raw(os.path.join(LAB,'best','new_instance_partial_39026.json'))
WIRE=set(u for u in range(L.NVARS) if v[u]==P)
print("variables currently equal to p:",len(WIRE))
# are they forced?  check the defining atom of each
byform=collections.Counter()
for u in sorted(WIRE):
    a=L.definer.get(u)
    byform['free' if a is None else 'gate'] += 1
print(" of which defined by a gate:",byform)
# reduce every atom mod p by killing monomials that touch a wire
dead=0; live=[]
for a in range(L.NA):
    keep={m:c for m,c in L.polys[a].items() if not (set(m)&WIRE)}
    keep={m:c%P for m,c in keep.items() if c%P}
    if not keep: dead+=1
    else: live.append((a,keep))
print(f"atoms identically 0 mod p (all monomials touch a p-wire): {dead} of {L.NA}")
print(f"atoms with real mod-p content: {len(live)}")
gl=[a for a,k in live if L.atom_out.get(a) is not None]
cl=[a for a,k in live if L.atom_out.get(a) is None]
print(f"   of them gates: {len(gl)}   checks: {len(cl)}")
# which variables survive mod p?
sv=set()
for a,k in live:
    for m in k: sv|=set(m)
print("variables appearing in the live mod-p system:",len(sv))
# how many live CHECK atoms are currently nonzero mod p at the checkpoint?
bad=[a for a,k in live if L.atom_out.get(a) is None and sum(c*__import__('math').prod([v[u] for u in m]) if m else c for m,c in k.items())%P]
print("live checks nonzero mod p at the checkpoint:",len(bad), bad[:20])
