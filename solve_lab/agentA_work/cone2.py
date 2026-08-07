"""Ancestor cone of a variable in the gate DAG, and whether a target var is inside."""
import sys, collections; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
def cone(t, maxn=200000):
    seen=set(); st=[t]
    while st:
        u=st.pop()
        if u in seen: continue
        seen.add(u)
        a=L.definer.get(u)
        if a is None: continue
        for w in L.avars[a]:
            if w!=u and w not in seen: st.append(w)
    return seen
for t in [int(x) for x in sys.argv[1:]]:
    c=cone(t)
    fr=[u for u in c if u not in L.definer]
    print('x%d cone size=%d free inputs=%d'%(t,len(c),len(fr)))
    print('   contains x14853=%s x14623=%s x9118=%s x8731=%s x24548=%s x25442=%s'%(
        14853 in c,14623 in c,9118 in c,8731 in c,24548 in c,25442 in c))
