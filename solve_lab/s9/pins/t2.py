import sys, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s9/pins')
from build import *
codes,_=H.load_equations()
d=pickle.load(open('atoms.pkl','rb')); src=d['atom_src']
atom2eq=pickle.load(open('atom2eq.pkl','rb'))
print(type(atom2eq))
for a in [21617,21619,21621,21623,29539,37662,37887,40826,19297,26733,36185,40812]:
    e=atom2eq.get(a,[]) if isinstance(atom2eq,dict) else []
    print(a, len(e), repr(src[a])[:200])
