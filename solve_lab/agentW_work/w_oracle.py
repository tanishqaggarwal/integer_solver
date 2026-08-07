"""Audit the consistency oracle O's negatives rest on.

solve_sparse returns None for FIVE distinct reasons, only three of which mean 'infeasible':
  'row X unsatisfiable' / 'row X: rhs %% c != 0' / 'core infeasible' / 'backsub non-integral'  -> genuinely infeasible
  'core too large' / 'coefficient blowup'                                                      -> GAVE UP
If any of O's 198,772 solves returned None for a give-up reason, that solve is not a negative.
"""
import sys, os, itertools, random, collections, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w_setup as S
sys.path.append('/home/user/integer_solver/solve_lab/agentE_work')
import sparse

def raw(Sset):
    return sparse.solve_sparse([S.rows[e] for e in Sset], [S.rhs[e] for e in Sset],
                               names=[str(x) for x in Sset], verbose=False, maxcore=600,
                               maxbits=10**7, maxcorebits=10**7)

msgs = collections.Counter()
t0=time.time(); n=0
# b=0 over every subset of FAIL (127), then a random sample of b=1 and b=2 across triples
for k in range(1,8):
    for P in itertools.combinations(S.FAIL,k):
        s,m,_ = raw(S.SAT+list(P)); n+=1
        msgs[m.split(':')[0].split(' unsat')[0] if s is None else 'OK'] += 1
random.seed(7)
trip = list(itertools.combinations(S.FAIL,3))
for _ in range(400):
    P = list(random.choice(trip))
    r1,r2 = random.sample(S.SAT,2)
    s,m,_ = raw([e for e in S.SAT if e not in (r1,r2)]+P); n+=1
    msgs[m.split(':')[0].split(' unsat')[0] if s is None else 'OK'] += 1
print('%d solves in %.1fs (%.1f/s)' % (n, time.time()-t0, n/(time.time()-t0)))
for k,v in msgs.most_common(): print('  %-40s %d' % (k,v))
