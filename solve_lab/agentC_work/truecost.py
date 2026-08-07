"""TRUE cost of a detach plan: every atom containing a detached handle variable becomes
nonzero (its value is a linear form in the detached values), so the cost is the union of
ALL their equations, minus the pins we satisfy by construction."""
import sys, json, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from close4 import *
W='/home/user/integer_solver/solve_lab/agentC_work/'
BI=json.load(open(W+'bitinfo.json'))
AE={a:frozenset(L.atom2eq.get(a,{})) for a in range(L.NA)}
eq2atoms=collections.defaultdict(set)
for a,es in AE.items():
    for e in es: eq2atoms[e].add(a)
rows=[]
for b,d in BI.items():
    b=int(b)
    H1=d['xpin']['H']; H2=d['ypin']['H']
    pins={d['xpin']['atom'],d['ypin']['atom']}
    S=set()
    for H in (H1,H2):
        S |= set(L.var_atoms[H])
    S-=pins
    # atoms that are gate definers of H itself are the detached ones (broken)
    E=frozenset().union(*[AE[a] for a in S]) if S else frozenset()
    # compensators: atoms fully inside E (extra free parameters if settable)
    cand=set()
    for e in E: cand|=eq2atoms[e]
    inside=[a for a in cand if AE[a]<=E]
    rows.append((len(E),len(S),len(inside),b,H1,H2,sorted(S)))
rows.sort()
print('TRUE cost table (|E_total| = union of equations of every atom touching the two handles):')
for r in rows[:20]:
    print('  |E|=%-3d natoms=%-3d inside=%-3d bit=x_%-6d H=(%d,%d)'%r[:6])
print()
print('distribution of |E_total|:',collections.Counter(r[0] for r in rows))
# calibrate on the deliverable: its broken atoms
D=[22229,22230,35758,35759,35760,35761,35762]
E=frozenset().union(*[AE[a] for a in D])
cand=set()
for e in E: cand|=eq2atoms[e]
inside=[a for a in cand if AE[a]<=E]
print('deliverable: |E|=%d  broken=%d  inside=%d  -> observed failing 7'%(len(E),len(D),len(inside)))
json.dump([[r[0],r[1],r[2],r[3],r[4],r[5],r[6]] for r in rows[:60]],open(W+'truecost.json','w'))
