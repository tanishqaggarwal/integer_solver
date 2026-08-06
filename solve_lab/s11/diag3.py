import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw, atomval
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
v = load_raw(os.path.join(LAB,'best','new_instance_partial_39026.json'))
def sh(x):
    s=str(x); return s if len(s)<26 else s[:12]+'..'+s[-8:]+f'<{len(s)}d>'
def nm(u):
    if v[u]==P: return f"x{u}=P"
    return f"x{u}={sh(v[u])}"
for a in list(range(22229,22236))+list(range(19087,19093))+list(range(10935,10939))+list(range(35756,35763)):
    Pp=L.polys[a]; out=L.atom_out.get(a)
    terms=' + '.join(f"{c}*{'*'.join('x%d'%u for u in m)}" for m,c in Pp.items())
    print(f"a{a:6d} out={str(out):14s} val={sh(atomval(a,v))}   {terms}")
    print(f"          vals: {[nm(u) for u in sorted(set(u for m in Pp for u in m))]}")
