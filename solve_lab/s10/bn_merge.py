"""bn_merge: merge sweep shards into a durable JSONL and analyse the cost law."""
import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'
sys.path.insert(0, HERE)
import bn_lib as B
import lib as L, tools as T, ad

bools=B.bools_map()
rows=[]
for i in range(4):
    p=os.path.join(HERE,f'bn_sweep_2_{i}.json')
    if os.path.exists(p): rows+=json.load(open(p))
rows.sort(reverse=True)
with open(os.path.join(HERE,'bn_sweep.jsonl'),'w') as f:
    for r in rows:
        f.write(json.dumps(dict(zip(['score','atom','var','boolc','oldval','nzatoms'],r)))+'\n')
print('merged',len(rows),'records -> bn_sweep.jsonl')
print('baseline 39026')
print('TOP 12:')
for r in rows[:12]:
    a=r[1]; print(f'   score={r[0]} a{a} x_{r[2]} boolc={r[3]} was={r[4]} nzatoms={r[5]} '
                  f'|E(a)|={len(L.atom2eq[a])}')
h=collections.Counter(r[0] for r in rows)
print('score histogram:',sorted(h.items(),reverse=True)[:12])
hn=collections.Counter(r[5] for r in rows)
print('nonzero-atom-count histogram:',sorted(hn.items())[:12])
# is cost exactly |E(a)| ?
ex=0; more=0
for r in rows:
    cost=39026-r[0]; ne=len(L.atom2eq[r[1]])
    if cost==ne: ex+=1
    elif cost>ne: more+=1
print(f'cost == |E(a)| for {ex}/{len(rows)};  cost > |E(a)| for {more}')
d=collections.Counter(39026-r[0]-len(L.atom2eq[r[1]]) for r in rows)
print('cost - |E(a)| histogram:',sorted(d.items())[:12])
# do the 7 baseline failures ever get fixed?
print()
best=rows[0]
print('best single flip:',best)
