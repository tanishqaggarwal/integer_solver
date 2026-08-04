import heal_harness as H
from collections import defaultdict, deque
p=H.p
d=H.loadd('best/new_instance_partial_39013.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward()

# For each equation, its free-input support (ancestors of its vars)
def eq_free_anc(i):
    s=set()
    for v in H.eqvars[i]:
        if v in H.freeinp: s.add(v)
        else: s |= H.anc.get(v,set())
    return s

eqfree=[None]*len(H.eqcode)
nnz=0
free2eq=defaultdict(list)
for i in range(len(H.eqcode)):
    s=eq_free_anc(i)
    eqfree[i]=s
    nnz+=len(s)
    for f in s: free2eq[f].append(i)
print("total incidence nnz (eq x free):", nnz)
print("free inputs appearing in >=1 eq:", len(free2eq))

# Connected components in bipartite graph; find comp of core frees
core_free={14853,12186,16742}
# also whatever free anc of x_3558 (second core) 
core_free |= H.anc.get(3558,set())
seen_eq=set(); seen_free=set()
dq=deque()
for f in core_free:
    if f in free2eq: seen_free.add(f); dq.append(('f',f))
while dq:
    typ,x=dq.popleft()
    if typ=='f':
        for i in free2eq[x]:
            if i not in seen_eq: seen_eq.add(i); dq.append(('e',i))
    else:
        for f in eqfree[x]:
            if f not in seen_free: seen_free.add(f); dq.append(('f',f))
print("coupling component: eqs=",len(seen_eq)," free=",len(seen_free))
# how many satisfied vs failing in component
F=set(H.fails())
print("failing eqs in component:", len(seen_eq&F), "of", len(F))
# nnz within component
compnnz=sum(len(eqfree[i]&seen_free) for i in seen_eq)
print("component nnz:", compnnz)
import pickle
pickle.dump({'eqfree':eqfree,'free2eq':dict(free2eq),'seen_eq':seen_eq,'seen_free':seen_free},
            open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/struct.pkl','wb'))
