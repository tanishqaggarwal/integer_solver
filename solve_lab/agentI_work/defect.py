import pickle, os, sys, collections
from model import Model
from fp import FpEngine, P
from boolscore import Fast
HERE = os.path.dirname(os.path.abspath(__file__))
d = pickle.load(open(os.path.join(HERE, 'bool_wit.pkl'), 'rb'))
val = d['val']; conf = d['conf']
F = Fast(); M = F.M; E = F.E
print("conflict atoms and their defect values:")
for a in conf:
    print(f"  a{a}: {M.src[a]}  ->  {E.eval_atom(a, [0 if x is None else x for x in val])}")
for v in [2287, 21889, 25156, 35389, 6671, 3023]:
    print(f"  X{v} = {val[v]}")
# cone of a17810
reason = d['val']  # placeholder
EOF = None
