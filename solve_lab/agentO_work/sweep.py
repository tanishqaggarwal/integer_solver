"""Agent O: representative sweep over the 106 proven pin-solvable bits, through E's
   LOG-16 simultaneous solve.  Read-only use of agentE_work modules."""
import sys, json, os, time, itertools, pickle, argparse
ED='/home/user/integer_solver/solve_lab/agentE_work'
sys.path.insert(0,ED)
os.chdir(ED)                      # E's modules load pickles by relative path
sys.set_int_max_str_digits(20_000_000)
import channels as C, engine as E
OD='/home/user/integer_solver/solve_lab/agentO_work'

FE=json.load(open(OD+'/feasbits.json'))
CH=json.load(open(ED+'/chan_cfg0.json'))['chan']
CHSET=[set(c) for c in CH]
FEAS=set(FE['A'])|set(FE['B'])
# feasible representatives per cfg0 channel
REPS=[sorted(CHSET[i]&FEAS) for i in range(3)]
INERT=sorted(FEAS-set().union(*CHSET))

def run(seed, maxr=3, maxv=2000):
    t0=time.time()
    try:
        r=C.simsolve(seed, maxr=maxr, maxv=maxv)
    except Exception as e:
        return None, f'ERR {type(e).__name__}: {e}', time.time()-t0
    if r is None: return None, 'nosol', time.time()-t0
    n,ns,av,v = r
    return (n,ns,av,v), 'ok', time.time()-t0
