import sys
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import harness as H
definer=H.definer; avars=H.avars; atoms=H.atoms; occ=H.occ
C=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
v=H.forward({18956:C})
seen=set()
def dump(u,depth,maxd):
    if u in seen or depth>maxd: return
    seen.add(u)
    dv=definer[u]
    print('  '*depth+f"x_{u} = {str(v[u])[:40]} :: {'FREE' if dv is None else atoms[dv[0]][:110]}")
    if dv is None: return
    for w in avars[dv[0]]:
        if w!=u: dump(w,depth+1,maxd)
for a in sys.argv[1:]:
    if a.startswith('d'): continue
maxd=int(sys.argv[-1])
for a in sys.argv[1:-1]:
    seen.clear(); dump(int(a),0,maxd)
    print()
