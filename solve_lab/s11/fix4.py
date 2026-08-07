import sys, os, json, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw, atomval
P=L.P
HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
v=load_raw(os.path.join(HERE,'data','fix2_round.json'))
v0=load_raw(os.path.join(LAB,'best','new_instance_partial_39026.json'))
def sh(x):
    s=str(x); return s if len(s)<22 else s[:9]+'..'+s[-5:]+f'<{len(s)}d>'
def trace(u, depth=0, seen=None):
    if seen is None: seen=set()
    if u in seen or depth>4: return
    seen.add(u)
    a=L.definer.get(u)
    pad='   '*depth
    if a is None:
        print(f"{pad}x{u} FREE = {sh(v[u])}  (was {sh(v0[u])})"); return
    Pp=L.polys[a]
    terms=' + '.join(f"{c}*{'*'.join('x%d'%t for t in m)}" for m,c in Pp.items())
    print(f"{pad}x{u} := via a{a}  val={sh(v[u])} (was {sh(v0[u])})  atomval={sh(atomval(a,v))}")
    print(f"{pad}    {terms[:200]}")
    if len(Pp)<8:
        for t in sorted(set(t for m in Pp for t in m)):
            if t!=u: trace(t, depth+1, seen)
for u in [14853,1308,24548,25442]:
    print("#########", u); trace(u); print()
