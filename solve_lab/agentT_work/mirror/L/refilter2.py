"""EXACT atom->equation incidence, built from the checker's own equations.

Each residual atom has exactly ONE free cofactor u (3,681 of them, verified in hcheck.py), and
u occurs nowhere else in the system.  So equation e contains atom a  <=>  u_a in vars(e).
That is exact and needs no model of the equation algebra.
"""
import sys, json, pickle, collections, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab')
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentT_work/mirror/F')
import checker as CK
from fwd import Engine,NV
from circ2 import vars_of
E=Engine()
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
defrhs={c[1]:c[2] for c in (E.cls[a] for a in E.order)}
Hd=pickle.load(open('handles.pkl','rb')); handle=set(Hd['handle'])
sys.setrecursionlimit(100000)
freeall={}
def fa(v):
    if v in freeall: return freeall[v]
    if v not in defrhs: freeall[v]={v}; return freeall[v]
    freeall[v]=set(); s=set()
    for u in vars_of(defrhs[v]): s|=fa(u)
    freeall[v]=s; return s
atomu={}
for a in E.res:
    s=set()
    for v in vars_of(E.atoms[a]): s|=fa(v)
    hs=[v for v in s if v in handle]
    if len(hs)==1: atomu[a]=hs[0]
print('atoms with a unique cofactor: %d'%len(atomu),flush=True)
t0=time.time(); CODES,VS=CK.load_equations(); print('checker eqs loaded %.0fs'%(time.time()-t0),flush=True)
u2eq=collections.defaultdict(set)
for i,vs in enumerate(VS):
    for v in vs:
        u2eq[v].add(i)
EQ={a:u2eq.get(u,set()) for a,u in atomu.items()}
DELIV=['((x7075*x8731)+x31864)','((5113045*(x7075*x9118))-x29854)',
       '((x4432-x19964)-x28730)','((x7068-x2099)-(7376877*x642))']
M25=set([2554,5324,6816,8124,8680,9041,9123,9421,11226,12231,12270,12350,14584,15558,18673,
         21000,22044,22534,22997,28929,29125,29330,32026,35512,38051])
u=set()
for a in DELIV: u|=EQ.get(a,set())
print('\nCALIBRATION with the EXACT map: deliverable 4 atoms touch %d equations'%len(u))
print('   vs M25: %d of 25 ; M25 not covered: %s'%(len(u&M25),sorted(M25-u)))
print('   deliverable checker failures subset? %s'%
      set([12231,12270,12350,14584,18673,22044,29125]).issubset(u))
pickle.dump({a:sorted(s) for a,s in EQ.items()},open('eqmap_exact.pkl','wb'))
# ---- my own baseline failing set, for the caveat M flagged ----
def exact_fail(vv):
    v=[0]*CK.NVARS; n=min(len(vv),CK.NVARS); v[:n]=vv[:n]
    ns={'v':v,'__builtins__':{}}
    return set(i for i,c in enumerate(CODES) if eval(c,ns)!=0)
D=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
vd=[0]*CK.NVARS
for k,x in D.items(): vd[int(k[2:])]=int(x)
# uncorrupted baseline = deliverable with its four corrupted handles (and their cofactors) zeroed
vb=list(vd)
for h in (642,28730,29854,31864,1329,9413,10903,17325,105,3387,5081,5676,11436,14393,14768,22820):
    vb[h]=0
B=exact_fail(vb)
print('\nMY baseline (deliverable with the 16 tuned handle/cofactor vars zeroed): %d failing'%len(B))
print('   |B & M25| = %d   B\\M25 = %s   M25\\B = %s'%(len(B&M25),sorted(B-M25)[:12],sorted(M25-B)[:12]))
USE=B&M25 if (B&M25) else M25
print('   FILTERING SET (intersection, per the caveat): %d equations'%len(USE))
json.dump(sorted(USE),open('filterset.json','w'))
