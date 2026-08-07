import sys, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import harness as H
definer=H.definer; avars=H.avars; atoms=H.atoms; occ=H.occ
C=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
v=H.forward({18956:C})
def show(u,depth=0,seen=None,maxd=4):
    if seen is None: seen=set()
    ind='  '*depth
    dv=definer[u]
    print(f"{ind}x_{u} = {str(v[u])[:50]} :: {'FREE' if dv is None else atoms[dv[0]][:110]}")
    print(f"{ind}   other atoms: {[ (i,atoms[i][:60]) for i in occ[u] if dv is None or i!=dv[0]][:6]}")
    if u in seen or depth>=maxd or dv is None: return
    seen.add(u)
    for w in avars[dv[0]]:
        if w!=u: show(w,depth+1,seen,maxd)
for aid in (20212,24403):
    print("="*100); print("ATOM",aid,atoms[aid])
    for u in avars[aid]: show(u,1,set(),2)
