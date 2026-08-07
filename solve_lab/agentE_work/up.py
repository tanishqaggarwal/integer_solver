import sys
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import harness as H
definer=H.definer; avars=H.avars; atoms=H.atoms; occ=H.occ
def chain(u,maxd=12,ind=0):
    seen=set()
    while True:
        dv=definer[u]
        print('  '*ind + f"x_{u} :: {'FREE' if dv is None else atoms[dv[0]][:100]}")
        if dv is None or u in seen: return
        seen.add(u)
        ups=[w for w in avars[dv[0]] if w!=u]
        if len(ups)==1: u=ups[0]; ind+=1; continue
        for w in ups: chain(w,maxd,ind+1)
        return
for u in [int(a) for a in sys.argv[1:]]:
    print("---",u); chain(u)
