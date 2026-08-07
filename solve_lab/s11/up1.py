import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
v=load_raw(os.path.join(LAB,'best','new_instance_partial_39026.json'))
def sh(x):
    s=str(x); return s if len(s)<20 else s[:8]+'..'+s[-4:]+f'<{len(s)}d>'
def trace(u, depth=0, seen=None, maxd=8):
    if seen is None: seen=set()
    pad='  '*depth
    a=L.definer.get(u)
    tag=f"x{u}({len(L.var_atoms[u])} atoms)"
    if a is None:
        print(f"{pad}{tag} FREE = {sh(v[u])}"); return
    if u in seen:
        print(f"{pad}{tag} ..."); return
    seen.add(u)
    Pp=L.polys[a]
    terms=' + '.join(f"{c}*{'*'.join('x%d'%t for t in m)}" for m,c in Pp.items())
    print(f"{pad}{tag} := a{a}  = {sh(v[u])}")
    print(f"{pad}   {terms[:160]}")
    if depth>=maxd or len(Pp)>6: return
    for t in sorted(set(t for m in Pp for t in m)):
        if t!=u: trace(t, depth+1, seen, maxd)
trace(2099)
