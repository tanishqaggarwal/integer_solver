import sys, os, json
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from ip7 import load_raw
from gmp1 import evalp, forwardp
P=L.P
sys.set_int_max_str_digits(400000)
LAB=os.path.join(HERE,'..')
def sh(x):
    s=str(x); return s if len(s)<14 else s[:6]+'..'+f'<{len(s)}d>'
for nm,f in [('checkpoint',os.path.join(LAB,'best','new_instance_partial_39026.json')),
             ('39018 (chan B)',os.path.join(HERE,'data','finish3_named.json'))]:
    v=[x%P for x in load_raw(f)]; forwardp(v)
    print(f"{nm}: x15298={sh(v[15298])} x5647={sh(v[5647])} x34606={sh(v[34606])} "
          f"x3896={sh(v[3896])} x38170={sh(v[38170])}")
    print(f"    x8599={sh(v[8599])} x21839={sh(v[21839])} x7304={sh(v[7304])} x25956={sh(v[25956])}")
for a,nm in [(23000,'OR-gate a23000')]:
    print(nm, L.polys[a])
# definitions
for u in [15298,5647,34606,3896,38170,8599,21839,7304,25956]:
    d=L.definer.get(u)
    t=(' + '.join(f"{c}*{'*'.join('x%d'%z for z in m)}" for m,c in L.polys[d].items())[:110]) if d is not None else 'FREE'
    print(f"  x{u} := {t}")
