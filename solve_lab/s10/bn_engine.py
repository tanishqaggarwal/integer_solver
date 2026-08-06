"""bn_engine: build the best boolean-carrier states and run the repair engine.

usage: bn_engine.py <var> <val>
"""
import os, sys, json, time
HERE='/home/user/integer_solver/solve_lab/s10'
sys.path.insert(0, HERE)
import bn_lib as B
import lib as L, tools as T, ad, engine

u=int(sys.argv[1]); x=int(sys.argv[2])
v0=L.load(B.BEST); BG=B.broken_gates(v0)
v=list(v0); v[u]=x
B.fwdb(v,BG,1)
s,fl,av=B.score(v)
print(f'carrier x_{u}={x}: start score {s}, failing {len(fl)}',flush=True)
T.save(v, os.path.join(HERE,f'bn_carrier_{u}_{x}.json'))

# engine with the broken gates protected: monkeypatch fwd to fwdb
_orig=ad.fwd
def fwdb_rounds(w, rounds=1):
    return B.fwdb(w, BG, 1)
ad.fwd=fwdb_rounds
engine.FORBID=set()
v2,cur=engine.run(v, f'bn_{u}_{x}', iters=40, budget=900)
ad.fwd=_orig
print('FINAL',cur[0],flush=True)
if cur[0]>39026:
    T.save(v2, os.path.join(HERE,'bn_best.json'))
    print('SAVED bn_best.json',flush=True)
