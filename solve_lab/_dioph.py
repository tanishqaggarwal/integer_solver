import heal_harness as H, json, re, random
from collections import defaultdict
p=H.p
VAR=re.compile(r'x_(\d+)')
d=H.loadd('best_agentA_39022.json')
base={v:d.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
F=[2554, 6816, 8124, 8680, 9421, 12231, 12270, 12350, 14584, 22044, 29125]
AFF=[1329, 1613, 6090, 7068, 8731, 9118, 10903, 17325, 21574, 27500]
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':H.val,'__builtins__':{}}
def setfree(v,x):
    H.val[v]=x
    for k in desc_of[v]: H.val[H.order[k]]=eval(H.gcode[k],ns)
def resid():
    return [eval(H.eqcode[i],ns) for i in F]
b0=resid()   # constant term b (residual at base)
# build A: A[k][j] = d resid_k / d AFF_j  (integer slope)
A=[[0]*len(AFF) for _ in F]
for j,v in enumerate(AFF):
    setfree(v, base[v]+1); r1=resid(); setfree(v, base[v])
    for k in range(len(F)): A[k][j]=r1[k]-b0[k]
# VERIFY joint affinity at a random integer point
random.seed(7)
pt={v: base[v]+random.randint(-5,5) for v in AFF}
for v in AFF: setfree(v, pt[v])
ract=resid()
for v in AFF: setfree(v, base[v])
ok=True
for k in range(len(F)):
    pred=b0[k]+sum(A[k][j]*(pt[AFF[j]]-base[AFF[j]]) for j in range(len(AFF)))
    if pred!=ract[k]: ok=False; print("NONLINEAR cross-term in eq",F[k],"pred!=act")
print("Joint affine model exact?",ok)
# Show the linear system A x = -b0 (we solve for delta = x - base)
print("\nResidual system (want A*delta = -b0), b0 bit-lengths:",[b.bit_length() for b in b0])
# Reduce mod p first: is A*delta ≡ -b0 mod p solvable? (necessary)
import sympy
Amp=sympy.Matrix([[A[k][j]%p for j in range(len(AFF))] for k in range(len(F))])
bmp=sympy.Matrix([(-b0[k])%p for k in range(len(F))])
# rank over GF(p)
Aug=Amp.row_join(bmp)
print("rank A mod p:", Amp.rank(iszerofunc=lambda x: x%p==0) if False else 'compute below')
json.dump({'A':[[A[k][j] for j in range(len(AFF))] for k in range(len(F))],
           'b0':[str(x) for x in b0],'AFF':AFF,'F':F,
           'base':{str(v):str(base[v]) for v in AFF}},
          open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/dioph.json','w'))
print("saved dioph.json")
