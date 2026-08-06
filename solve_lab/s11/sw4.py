import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
v=load_raw(os.path.join(HERE,'data','sw3_out.json'))
def sh(x):
    s=str(x); return s if len(s)<20 else s[:8]+'..'+s[-4:]+f'<{len(s)}d>'
def trace(u,d=0,seen=None,maxd=6):
    seen=seen if seen is not None else set()
    pad='  '*d; a=L.definer.get(u)
    if a is None:
        print(f"{pad}x{u} FREE={sh(v[u])} atoms={len(L.var_atoms[u])} %p={'0' if v[u]%P==0 else 'nz'}"); return
    if u in seen: print(f"{pad}x{u} ..."); return
    seen.add(u)
    Pp=L.polys[a]
    print(f"{pad}x{u}:=a{a} ={sh(v[u])} %p={'0' if v[u]%P==0 else 'nz'}  "+
          ' + '.join(f"{c}*{'*'.join('x%d'%t for t in m)}" for m,c in Pp.items())[:130])
    if d>=maxd or len(Pp)>5: return
    for t in sorted(set(t for m in Pp for t in m)):
        if t!=u: trace(t,d+1,seen,maxd)
for u in [2239,31731,9106]:
    print("#####",u); trace(u); print()
