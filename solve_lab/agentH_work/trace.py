import ev, pickle, sys
from collections import deque
F=ev.F; definer=F['definer']; atom_src=ev.atom_src; atom_vars=ev.atom_vars
def cone(v, maxn=200):
    seen=set(); Q=deque([v]); out=[]
    while Q:
        u=Q.popleft()
        if u in seen: continue
        seen.add(u)
        a=definer[u]
        if a<0:
            out.append((u,None)); continue
        out.append((u,atom_src[a]))
        for w in atom_vars[a]:
            if w!=u and w not in seen: Q.append(w)
        if len(seen)>maxn: break
    return out
for v in [int(x) for x in sys.argv[1:]]:
    print('=== cone of x_%d ==='%v)
    for u,s in cone(v):
        print('  x_%d ='%u, (s[:220] if s else 'FREE INPUT'))
