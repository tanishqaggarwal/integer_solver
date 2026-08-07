"""S11 step 56: can x_14623 be DECOUPLED from x_33462 inside the K1 web?

a21617 needs x_14623 to move mod p; the K1 web ties it to x_33462, which a31672
pins to K1.  A uniform shift of the whole class fixes a21617 but breaks a31672
(net 0).  A PARTIAL shift needs only the gadgets linking the two -- price them.
"""
import os, sys, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, ad
P = ad.P
v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
K1 = 97171863764434070215824145711260403004952728652948669662983319257693684265837195009100680
cls = [u for u in range(L.NVARS) if v[u] == K1]
print(f'K1 class: {len(cls)} variables; free inputs among them: '
      f'{[u for u in cls if u not in L.definer]}')
# which atoms link two members of the class?
links = collections.defaultdict(list)
for a in range(L.NA):
    mem = [u for u in set(L.avars[a]) if v[u] == K1]
    if len(mem) >= 2: links[tuple(sorted(mem))].append(a)
print(f'atoms containing >=2 class members: {sum(len(x) for x in links.values())}')
# build the link graph and find the path from x_14623 to x_33462
adj = collections.defaultdict(set)
edge_atom = {}
for pair, ats in links.items():
    for i in range(len(pair)):
        for j in range(i+1, len(pair)):
            adj[pair[i]].add(pair[j]); adj[pair[j]].add(pair[i])
            edge_atom[(pair[i], pair[j])] = ats
            edge_atom[(pair[j], pair[i])] = ats
import collections as C
q = C.deque([(14623, [14623])]); seen = {14623}; path = None
while q:
    u, pth = q.popleft()
    if u == 33462: path = pth; break
    for w in adj[u]:
        if w not in seen: seen.add(w); q.append((w, pth + [w]))
print(f'\npath from x_14623 to x_33462 in the link graph: {path}')
if path:
    tot = 0
    for i in range(len(path) - 1):
        ats = edge_atom[(path[i], path[i+1])]
        cost = min(len(L.atom2eq[a]) for a in ats)
        cheap = [a for a in ats if len(L.atom2eq[a]) == cost]
        tot += cost
        print(f'  x_{path[i]} -- x_{path[i+1]}  via {[(a, len(L.atom2eq[a])) for a in ats]}'
              f'   cheapest {cost}')
    print(f'  cheapest single edge to cut on this path: '
          f'{min(min(len(L.atom2eq[a]) for a in edge_atom[(path[i], path[i+1])]) for i in range(len(path)-1))}')
print(f'\nreachable class members from x_14623: {len(seen)} of {len(cls)}')
print(f'  is x_33462 reachable? {33462 in seen}')
