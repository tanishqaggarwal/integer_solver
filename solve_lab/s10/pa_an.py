"""Analyse the whole-instance break census: which atoms can be broken cheaply, and
which reachable atom supports S give the best balance-law score |E(S)|-|S|+c."""
import os, sys, collections, json
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L
d={}
for f in ('pa_census_a.json','pa_census_b.json'):
    d.update(json.load(open(os.path.join(HERE,f))))
print('census entries',len(d))
BASE={22229,22230,35758,35759,35760,35761,35762}
h=collections.Counter(r['fail'] for r in d.values())
print('failing-count histogram:',sorted(h.items())[:20])
# supports reached
sup=collections.Counter()
for t,r in d.items():
    sup[tuple(r['nz'])]+=1
print('distinct supports reached:',len(sup))
# best NEW supports (different from BASE)
rows=[]
for t,r in d.items():
    S=set(r['nz'])
    if S==BASE: continue
    E=set()
    for a in S: E|=set(L.atom2eq[a])
    rows.append((r['fail'],len(E)-len(S),len(S),len(E),int(t),r['a'],tuple(sorted(S-BASE)),tuple(sorted(BASE-S))))
rows.sort()
print('\nBEST measured relocations (fail = actual failing count):')
seen=set()
for fail,defi,ns,ne,t,a,extra,gone in rows:
    k=(extra,gone)
    if k in seen: continue
    seen.add(k)
    print(f'  fail={fail:<4} |S|={ns} |E|={ne} deficiency={defi}  x_{t} (a{a})  +{list(extra)[:8]} -{list(gone)}')
    if len(seen)>=30: break
# minimum deficiency achieved over all reached supports
rows.sort(key=lambda z:z[1])
print('\nLOWEST deficiency |E(S)|-|S| over all reached supports:')
seen=set()
for fail,defi,ns,ne,t,a,extra,gone in rows[:4000]:
    k=(extra,gone)
    if k in seen: continue
    seen.add(k)
    print(f'  deficiency={defi} fail={fail} |S|={ns} |E|={ne}  x_{t} (a{a}) +{list(extra)[:8]} -{list(gone)}')
    if len(seen)>=20: break
# cheapest single extra atom broken
best={}
for fail,defi,ns,ne,t,a,extra,gone in rows:
    if len(extra)==1 and not gone:
        b=extra[0]
        if b not in best or fail<best[b][0]: best[b]=(fail,t,ne)
print('\ncheapest single-atom additions (cost above 7):')
for b,(fail,t,ne) in sorted(best.items(),key=lambda kv:kv[1][0])[:25]:
    print(f'  a{b} fp={len(L.atom2eq[b])} -> total failing {fail} (delta {fail-7}) via x_{t}')
