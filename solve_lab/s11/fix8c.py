import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
v=load_raw(os.path.join(LAB,'best','new_instance_partial_39026.json'))
def sh(x):
    s=str(x); return s if len(s)<20 else s[:8]+'..'+s[-4:]+f'<{len(s)}d>'
for a in [3575,3576,36602]:
    Pp=L.polys[a]
    print(f"a{a} out={L.atom_out.get(a)} nterms={len(Pp)} neq={len(L.atom2eq.get(a,{}))} val={sh(L.evalpoly(Pp,v))}")
    if len(Pp)<10:
        for m,c in Pp.items(): print("   ",c,'*','*'.join('x%d'%u for u in m),[sh(v[u]) for u in m])
    print()
# trace x26777 upstream
def trace(u, d=0, seen=None):
    seen=seen or set()
    if u in seen or d>6: return
    seen.add(u)
    a=L.definer.get(u); pad='  '*d
    if a is None: print(f"{pad}x{u} FREE={sh(v[u])} ({len(L.var_atoms[u])} atoms)"); return
    Pp=L.polys[a]
    print(f"{pad}x{u}:=a{a} ({len(L.var_atoms[u])} atoms) ={sh(v[u])}  "+
          ' + '.join(f"{c}*{'*'.join('x%d'%t for t in m)}" for m,c in Pp.items())[:150])
    if len(Pp)>5: return
    for t in sorted(set(t for m in Pp for t in m)):
        if t!=u: trace(t,d+1,seen)
trace(26777)
