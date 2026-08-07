import sys, pickle, json, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import harness as H
C=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
base={18956:C}
def cone(u0):
    seen=set(); st=[u0]
    while st:
        u=st.pop()
        if u in seen: continue
        seen.add(u)
        dv=H.definer[u]
        if dv is None: continue
        for w in H.avars[dv[0]]:
            if w!=u: st.append(w)
    return seen
res={}
for root in (7715,34554):
    cn=cone(root)
    frees=[u for u in cn if H.definer[u] is None]
    print("root",root,"free in cone",len(frees))
    for f in frees:
        s=dict(base); s[f]=1
        v=H.forward(s); ff,av=H.eqfails(v)
        res[(root,f)]=(v[7715],v[34554],len(ff),sorted(av))
        if v[root]==1:
            print(f"  x_{f}=1 -> x_7715={v[7715]} x_34554={v[34554]} fails={len(ff)} bad={sorted(av)}")
pickle.dump(res,open('sel.pkl','wb'))
best=sorted(res.items(), key=lambda kv: kv[1][2])[:15]
print("--- best by fails ---")
for k,val in best: print(k,val)
