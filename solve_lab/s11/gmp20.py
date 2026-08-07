import sys, os, json, pickle
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from gmp1 import evalp, forwardp
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
D=pickle.load(open(os.path.join(HERE,'data','resp_modp.pkl'),'rb'))
cols=D['cols']; base=D['base']; bd=D['bd']
def sh(x):
    s=str(x); return s if len(s)<16 else s[:7]+'..'+f'<{len(s)}d>'
for a in [3568,3578,36040,29253]:
    movers=[u for u,d in cols.items() if a in d]
    Pp=L.polys[a]
    print(f"a{a}: moved by {len(movers)} knobs {movers[:8]}")
    for m,c in Pp.items():
        print(f"     {sh(c)} * {'*'.join('x%d'%u for u in m)}  vals={[sh(base[u]) for u in m]}")
    # the handle
    hs=[m[0] for m,c in Pp.items() if len(m)==1 and abs(c)<10**12]
    for h in hs:
        d=L.definer.get(h)
        t=(' + '.join(f"{cc}*{'*'.join('x%d'%z for z in mm)}" for mm,cc in L.polys[d].items())[:90]) if d is not None else 'FREE'
        print(f"     handle x{h}: {t}")
