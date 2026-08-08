import sys, os, time
sys.path.insert(0, os.path.abspath('../..'))
import numpy as np
from synth.solver.model import build_modmul
from run_modmul import enumerate_ground_states
from sdp import qubo_to_C, sdp_min, rank_of
s=int(sys.argv[1])
t=time.time(); model=build_modmul(s,mode='wallace',seed=3); Q=model['Q']
states,pinned=enumerate_ground_states(model)
t1=time.time()-t
t=time.time(); C=qubo_to_C(Q.Q,Q.n); val,X,V=sdp_min(C,seed=1); t2=time.time()-t
print(f"s={s} n={Q.n} ngs={len(states)} pinned={len(pinned)} SDPopt={val:+.4f} enum_t={t1:.1f}s sdp_t={t2:.1f}s")
