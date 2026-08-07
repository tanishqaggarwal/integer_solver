import sys, pickle, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s9/eff')
import lib as L
W='/home/user/integer_solver/solve_lab/agentC_work/'
D=pickle.load(open(W+'supp3.pkl','rb'))
LN=pickle.load(open(W+'lin.pkl','rb'))['lin']
outs=D['outs']; free=D['free']; topo=D['topo']; csupp=D['csupp']; supp=D['supp']
freeset=set(free)
P=2**256-2**32-977
CHK=[a for a in range(L.NA) if a not in L.atom_out]
# order of gate evaluation restricted to defined vars
GORDER=[u for u in topo if u in LN]
def forward(v):
    for u in GORDER:
        c,rest=LN[u]; s=0
        for cc,m in rest:
            t=cc
            for w in m: t*=v[w]
            s+=t
        v[u]= -s//c if s%c==0 else 0
    return v
def evalatom(a,v):
    s=0
    for m,c in L.polys[a].items():
        t=c
        for w in m: t*=v[w]
        s+=t
    return s
def score(v):
    av=L.all_atom_values(v)
    return L.NEQ-len(L.failing_eqs(av)), av
