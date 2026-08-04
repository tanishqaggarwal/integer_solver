import heal_harness as H, re, json, pickle, time
p=H.p
AD=pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/adcache.pkl','rb'))
gate_parsed=AD['gate_parsed']; val=AD['val']
for v in range(H.NVARS): H.val[v]=val[v]
# forward-mode tangent: given seed free f (dv[f]=1), propagate dv through gates
def tangent(f):
    dv=[0]*H.NVARS
    dv[f]=1
    for t in H.order:
        pr=gate_parsed[t]; op=pr[0]
        if op=='id':
            k=pr[1]; dv[t]=dv[k[1]] if k[0]=='v' else 0
        elif op=='add':
            k1,k2=pr[1],pr[2]
            dv[t]=((dv[k1[1]] if k1[0]=='v' else 0)+(dv[k2[1]] if k2[0]=='v' else 0))%p
        elif op=='sub':
            k1,k2=pr[1],pr[2]
            dv[t]=((dv[k1[1]] if k1[0]=='v' else 0)-(dv[k2[1]] if k2[0]=='v' else 0))%p
        elif op=='mul':
            k1,k2=pr[1],pr[2]
            v1=val[k1[1]] if k1[0]=='v' else k1[1]
            v2=val[k2[1]] if k2[0]=='v' else k2[1]
            d1=dv[k1[1]] if k1[0]=='v' else 0
            d2=dv[k2[1]] if k2[0]=='v' else 0
            dv[t]=(v1*d2+v2*d1)%p
    return dv
# eq derivative via atom structure: load atoms and eq->atom map with coeffs?
# We don't have coeffs. Instead eval eq residual derivative by dual: use eqvars to compute
# d_eq = sum_{v in eqvars} (partial eq/partial v)*dv[v].
# partial eq/partial v: compute via central-ish exact? eq is polynomial; use dual per eqvar is slow.
# Better: eval eq with all-vars dual using the tangent dv directly (one pass) via compiled code with dual v.
# Implement Dual and eval closure eqs.
class D:
    __slots__=('a','b')
    def __init__(s,a,b=0): s.a=a%p; s.b=b%p
    def __add__(s,o):
        o=o if isinstance(o,D) else D(o)
        return D(s.a+o.a, s.b+o.b)
    __radd__=__add__
    def __sub__(s,o):
        o=o if isinstance(o,D) else D(o); return D(s.a-o.a, s.b-o.b)
    def __rsub__(s,o):
        o=o if isinstance(o,D) else D(o); return D(o.a-s.a, o.b-s.b)
    def __mul__(s,o):
        o=o if isinstance(o,D) else D(o); return D(s.a*o.a, s.a*o.b+s.b*o.a)
    __rmul__=__mul__
    def __neg__(s): return D(-s.a,-s.b)
# baseline fails / satisfied
Fset=set(H.fails())
sat=[i for i in range(len(H.eqcode)) if i not in Fset]
print(f"satisfied {len(sat)}, failing {len(Fset)}")
# closure BFS: start frees = movers of S,T
gS=AD['gS']; gT=AD['gT']
R=set(f for f in H.freeinp if gS[f]%p or gT[f]%p)
print(f"initial movers {len(R)}")
# For efficiency, build var->eqs incidence for satisfied eqs
from collections import defaultdict
var2eq=defaultdict(list)
for i in sat:
    for v in H.eqvars[i]:
        var2eq[v].append(i)
pickle.dump({'sat':sat,'var2eq':dict(var2eq)}, open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/satcache.pkl','wb'))
# We'll compute the Jacobian columns for frees in R, over guard eqs that are fragile.
# First find fragile guards: perturb all R together with random tangent -> which sat eqs get nonzero deriv.
import random
random.seed(1)
# combined tangent = sum of random*tangent(f) is expensive (len(R) passes). Do it.
t0=time.time()
combo=[0]*H.NVARS
coefs={f:random.randrange(1,p) for f in R}
# propagate combined seed at once: set dv[f]=coefs[f], propagate
dv=[0]*H.NVARS
for f in R: dv[f]=coefs[f]
for t in H.order:
    pr=gate_parsed[t]; op=pr[0]
    if op=='id':
        k=pr[1]; dv[t]=dv[k[1]] if k[0]=='v' else 0
    elif op=='add':
        k1,k2=pr[1],pr[2]; dv[t]=((dv[k1[1]] if k1[0]=='v' else 0)+(dv[k2[1]] if k2[0]=='v' else 0))%p
    elif op=='sub':
        k1,k2=pr[1],pr[2]; dv[t]=((dv[k1[1]] if k1[0]=='v' else 0)-(dv[k2[1]] if k2[0]=='v' else 0))%p
    elif op=='mul':
        k1,k2=pr[1],pr[2]
        v1=val[k1[1]] if k1[0]=='v' else k1[1]; v2=val[k2[1]] if k2[0]=='v' else k2[1]
        d1=dv[k1[1]] if k1[0]=='v' else 0; d2=dv[k2[1]] if k2[0]=='v' else 0
        dv[t]=(v1*d2+v2*d1)%p
# now eval each sat eq derivative using dual: v[i]=D(val[i],dv[i]); eq deriv = .b
vdual=[D(val[i],dv[i]) for i in range(H.NVARS)]
ns={'v':vdual,'__builtins__':{}}
fragile=[]
for i in sat:
    r=eval(H.eqcode[i],ns)
    if isinstance(r,D):
        if r.b%p!=0: fragile.append(i)
    # if not D, derivative 0
print(f"fragile guards under random R-combo: {len(fragile)} (time {time.time()-t0:.1f}s)")
print("sample fragile:", sorted(fragile)[:20])
pickle.dump({'R':sorted(R),'fragile':sorted(fragile)}, open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/fragcache.pkl','wb'))
