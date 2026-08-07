"""For each candidate broken-atom seed set, compute E = union of its equations and
n = number of atoms whose entire equation set lies inside E (free compensators).
balance law:  failing = |E| - n + c."""
import sys, json, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from close4 import *
W='/home/user/integer_solver/solve_lab/agentC_work/'
BI=json.load(open(W+'bitinfo.json'))
AE={a:frozenset(L.atom2eq.get(a,{})) for a in range(L.NA)}
eq2atoms=collections.defaultdict(set)
for a,es in AE.items():
    for e in es: eq2atoms[e].add(a)
def inside(E):
    cand=set()
    for e in E: cand|=eq2atoms[e]
    return [a for a in cand if AE[a]<=E]
rows=[]
seen=set()
for b,d in BI.items():
    for seed in [(d['hx'],),(d['hy'],),(d['hx'],d['hy'])]:
        if seed in seen: continue
        seen.add(seed)
        E=frozenset().union(*[AE[a] for a in seed])
        S=inside(E)
        rows.append((len(E)-len(S),len(E),len(S),seed,int(b)))
rows.sort()
print('best 25 clusters by |E| - n :')
for r in rows[:25]:
    print('   |E|-n=%-4d |E|=%-3d n=%-3d seed=%s bit=x_%d'%r)
# also the deliverable's cluster
D=(22229,22230,35758,35759,35760,35761,35762)
E=frozenset().union(*[AE[a] for a in D]); S=inside(E)
print('deliverable cluster: |E|=%d n=%d  |E|-n=%d  atoms=%s'%(len(E),len(S),len(E)-len(S),sorted(S)))
# global scan: every atom as a seed
rows2=[]
for a in range(L.NA):
    E=AE[a]
    if not E or len(E)>20: continue
    S=inside(E)
    rows2.append((len(E)-len(S),len(E),len(S),a))
rows2.sort()
print('global single-atom seeds, best 20:')
for r in rows2[:20]:
    print('   |E|-n=%-4d |E|=%-3d n=%-3d a%d  %s'%(r[0],r[1],r[2],r[3],L.atom_src[r[3]][:80]))
json.dump([[r[0],r[1],r[2],r[3]] for r in rows2[:400]],open(W+'clusters.json','w'))
