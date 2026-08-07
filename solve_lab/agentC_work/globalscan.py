"""CORRECTED objective.  A cluster is good iff |E(S)| - (#INDEPENDENTLY SETTABLE atoms in S) is
small.  Shadow atoms (dependent linear combinations) inflate |S| without adding freedom, so the
right count is the number of handle-definition atoms (each carries its own free p-multiple knob).
Deliverable: |E|=12, settable=7, c=2 -> failing 7.  Beat it => |E| - settable + c <= 6."""
import sys, json, collections, fractions
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from close4 import *
W='/home/user/integer_solver/solve_lab/agentC_work/'
AE={a:frozenset(L.atom2eq.get(a,{})) for a in range(L.NA)}
eq2atoms=collections.defaultdict(set)
for a,es in AE.items():
    for e in es: eq2atoms[e].add(a)
v0=[0]*L.NVARS
_o2,_F,_fw,_C=mk(())
_fw(v0)
WIR={u for u in range(L.NVARS) if v0[u]==P}
FREE0=set(u for u in range(L.NVARS) if u not in outs)
# handle-definition atoms: out = wire * free   (settable: detach out, vary the free var by 1 -> value moves by p)
HD=set()
for a,(c,t) in L.atom_out.items():
    Pp=L.polys[a]
    if len(Pp)!=2: continue
    q=[m for m in Pp if len(m)==2]
    if len(q)!=1: continue
    w=set(q[0])
    if len(w&WIR)==1 and len((w-WIR)&FREE0)==1: HD.add(a)
print('settable handle-definition atoms:',len(HD))
def rk(mat):
    A=[[fractions.Fraction(x) for x in r] for r in mat]
    n=len(A); m=len(A[0]) if n else 0; r=0
    for c in range(m):
        piv=None
        for i in range(r,n):
            if A[i][c]!=0: piv=i; break
        if piv is None: continue
        A[r],A[piv]=A[piv],A[r]; pv=A[r][c]; A[r]=[x/pv for x in A[r]]
        for i in range(n):
            if i!=r and A[i][c]!=0:
                f=A[i][c]; A[i]=[x-f*y for x,y in zip(A[i],A[r])]
        r+=1
    return r
def score_cluster(seed):
    E=frozenset().union(*[AE[a] for a in seed])
    cand=set()
    for e in E: cand|=eq2atoms[e]
    inside=[a for a in cand if AE[a]<=E]
    sett=[a for a in inside if a in HD]
    if not sett: return None
    M=[[L.eq_atoms[e][2].get(a,0) for a in sett] for e in sorted(E)]
    r=rk(M)
    # max equations satisfiable <= rank of the settable block's column space seen by E, minus congruences
    return (len(E)-r, len(E), len(sett), r, tuple(sorted(seed)))
D=[22229,22230,35758,35759,35760,35761,35762]
print('deliverable cluster:',score_cluster(D),' (observed failing 7)')
rows=[]
for a in HD:
    s=score_cluster((a,))
    if s and s[1]<=20: rows.append(s)
rows.sort()
print('best 25 single-seed clusters by |E| - rank(settable block):')
for r in rows[:25]:
    print('   |E|-rank=%-3d |E|=%-3d settable=%-3d rank=%-3d seed=%s'%r)
json.dump([list(r[:4])+[list(r[4])] for r in rows[:300]],open(W+'globalscan.json','w'))
