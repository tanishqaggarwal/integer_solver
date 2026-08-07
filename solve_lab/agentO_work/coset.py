"""The deliverable's residual lives in exactly 12 equations over 8 atoms.  Enumerate every
   subset of those 12 that a NONZERO residual vector can satisfy, to bound the achievable
   failing-equation count from below."""
import sys, json, collections, itertools
from fractions import Fraction
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentO_work')
import simO, engine as E, harness as H
d=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
vd=[0]*E.NV
for k,x in d.items(): vd[int(k.split('_')[1])]=int(x)
bad=E.badatoms(vd); BAD=sorted(bad)
touch=collections.defaultdict(dict)
for e,(issq,outer,terms) in enumerate(H.eqt):
    for c,a in terms:
        if a in bad: touch[e][a]=c
EQ=sorted(touch)
M=[[touch[e].get(a,0) for a in BAD] for e in EQ]
print('atoms',BAD); print('eqs',EQ)
# identical columns?
for i in range(len(BAD)):
    for j in range(i+1,len(BAD)):
        if all(r[i]==r[j] for r in M): print(f'  columns for atoms {BAD[i]} and {BAD[j]} are IDENTICAL')

def rank(rows):
    R=[[Fraction(x) for x in r] for r in rows]; n=len(R); m=len(R[0]) if n else 0
    rk=0
    for c in range(m):
        p=None
        for i in range(rk,n):
            if R[i][c]: p=i; break
        if p is None: continue
        R[rk],R[p]=R[p],R[rk]
        pv=R[rk][c]
        for i in range(n):
            if i!=rk and R[i][c]:
                f=R[i][c]/pv
                for k in range(c,m): R[i][k]-=f*R[rk][k]
        rk+=1
    return rk
full=rank(M); print('rank of full 12x8:',full,'-> only r=0 satisfies all 12' if full==len(BAD) else '')
NC=len(BAD)
best=[]
for mask in range(1<<len(EQ)):
    S=[i for i in range(len(EQ)) if mask>>i&1]
    if len(S)<=len(best): continue
    if rank([M[i] for i in S])<NC:   # nonzero solution exists over Q
        best=S
print('max #equations satisfiable by a NONZERO residual on these 8 atoms:',len(best))
print('  that set:',[EQ[i] for i in best],'-> failing >=',len(EQ)-len(best))
# all maximal sets of that size
allmax=[]
for S in itertools.combinations(range(len(EQ)),len(best)):
    if rank([M[i] for i in S])<NC: allmax.append([EQ[i] for i in S])
print('  count of such sets:',len(allmax))
for s in allmax[:12]: print('   ',s)
cur=[e for e in EQ if sum(c*bad[a] for a,c in touch[e].items())==0]
print('currently satisfied:',cur,'(%d)'%len(cur))
