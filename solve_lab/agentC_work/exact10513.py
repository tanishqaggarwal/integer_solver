"""EXACT endgame on the x_10513 cluster.
1. build the detach state; 2. identify cluster E (equations) and S (atoms inside E);
3. MEASURE the achievable change-lattice on S by probing every free input (delta 1 and delta p),
   keeping only moves whose atom-change support lies inside S;
4. enumerate subsets of E and test exact integer solvability -> minimum failing count."""
import sys, json, time, collections, itertools
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from close4 import *
W='/home/user/integer_solver/solve_lab/agentC_work/'
TC={r[3]:r for r in json.load(open(W+'truecost.json'))}
b=10513
row=TC[b]; H1,H2=row[4],row[5]
AE={a:frozenset(L.atom2eq.get(a,{})) for a in range(L.NA)}
eq2atoms=collections.defaultdict(set)
for a,es in AE.items():
    for e in es: eq2atoms[e].add(a)
Sseed=set(row[6])
E=frozenset().union(*[AE[a] for a in Sseed])
cand=set()
for e in E: cand|=eq2atoms[e]
S=sorted(a for a in cand if AE[a]<=E)
E=sorted(E)
print('cluster: |E|=%d  |S|=%d'%(len(E),len(S)))
for a in S: print('   a%-6d eqs=%-2d  %s'%(a,len(AE[a]),L.atom_src[a][:110]))
M=[[L.eq_atoms[e][2].get(a,0) for a in S] for e in E]
print('M (equation x atom):')
for e,r in zip(E,M): print('   eq%-6d sq=%-5s %s'%(e,L.eq_atoms[e][1],r))
# rank over Q
def rank(mat):
    import fractions
    A=[[fractions.Fraction(x) for x in r] for r in mat]
    n=len(A); m=len(A[0]) if n else 0; r=0
    for c in range(m):
        piv=None
        for i in range(r,n):
            if A[i][c]!=0: piv=i; break
        if piv is None: continue
        A[r],A[piv]=A[piv],A[r]
        pv=A[r][c]
        A[r]=[x/pv for x in A[r]]
        for i in range(n):
            if i!=r and A[i][c]!=0:
                f=A[i][c]; A[i]=[x-f*y for x,y in zip(A[i],A[r])]
        r+=1
    return r
print('rank(M) over Q =',rank(M))
