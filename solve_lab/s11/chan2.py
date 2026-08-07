import sys, os, json
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from ip7 import load_raw
from gmp1 import evalp, forwardp
P=L.P; sys.set_int_max_str_digits(400000); LAB=os.path.join(HERE,'..')
v=[x%P for x in load_raw(os.path.join(LAB,'best','new_instance_partial_39026.json'))]; forwardp(v)
def sh(x):
    s=str(x); return s if len(s)<12 else s[:5]+'..'
def trace(u,d=0,seen=None):
    seen=seen if seen is not None else set()
    if u in seen or d>7: return
    seen.add(u); pad='  '*d
    a=L.definer.get(u)
    if a is None: print(f"{pad}x{u} FREE = {sh(v[u])}"); return
    Pp=L.polys[a]
    print(f"{pad}x{u}={sh(v[u])} := "+' + '.join(f"{c}*{'*'.join('x%d'%t for t in m)}" for m,c in Pp.items())[:100])
    if len(Pp)>4: return
    for t in sorted(set(t for m in Pp for t in m)):
        if t!=u: trace(t,d+1,seen)
for u in [7715,34554,23597,19271]:
    print("=====",u); trace(u)
