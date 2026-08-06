import heal_harness as H, json, re
from collections import defaultdict
p=H.p
VAR=re.compile(r'x_(\d+)')
d=H.loadd('best_agentA_39022.json')
base={v:d.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
F=[2554, 6816, 8124, 8680, 9421, 12231, 12270, 12350, 14584, 22044, 29125]
RF=[950, 1329, 1613, 3629, 4432, 6090, 6947, 7068, 8731, 8976, 9118, 9413, 10422, 10903, 11099, 15120, 15324, 17325, 21574, 22526, 27500, 33168, 34868, 35531]
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':H.val,'__builtins__':{}}
def setfree(v,x):
    H.val[v]=x
    for k in desc_of[v]: H.val[H.order[k]]=eval(H.gcode[k],ns)
def resid_modp():
    return [eval(H.eqcode[i],ns)%p for i in F]
b0=resid_modp()
# full Jacobian mod p (central-ish difference, but mod p use +1 slope from analytic? use +1 and +2 to get slope; but nonlinear -> slope at base via (f(x+1)-f(x)) is not exact deriv. Use symbolic-free: derivative mod p via f(x+1)-f(x) minus curvature. Better: use two-point to get 1st-order + treat. For Newton mod p we want the analytic Jacobian.)
# Analytic Jacobian: J_kj = (f(base+e_j)-f(base-e_j))/2 mod p won't be exact either for nonlinear.
# Instead compute exact partial derivative via: perturb by t, fit? For a low-degree poly, deriv at base = coefficient. Use 3 points to extract linear coeff exactly:
# f(x0+h) = f0 + h*f' + h^2/2 f'' ...; use f(x0+1)-f(x0-1) = 2 f' + (2/6)f''' ... For quadratic exact: (f(+1)-f(-1))/2 = f'. 
J=[[0]*len(RF) for _ in F]
for j,v in enumerate(RF):
    setfree(v, base[v]+1); rp=resid_modp()
    setfree(v, base[v]-1); rm=resid_modp()
    setfree(v, base[v])
    inv2=pow(2,p-2,p)
    for k in range(len(F)):
        J[k][j]=((rp[k]-rm[k])*inv2)%p
# rank of J and [J|-b0] mod p
def rank_modp(M):
    M=[row[:] for row in M]; rows=len(M); cols=len(M[0]) if rows else 0; r=0; pivcols=[]
    for c in range(cols):
        piv=None
        for i in range(r,rows):
            if M[i][c]%p!=0: piv=i;break
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]; inv=pow(M[r][c]%p,p-2,p)
        M[r]=[(x*inv)%p for x in M[r]]
        for i in range(rows):
            if i!=r and M[i][c]%p!=0:
                f=M[i][c]; M[i]=[(M[i][t]-f*M[r][t])%p for t in range(cols)]
        pivcols.append(c); r+=1
    return r,pivcols
rJ,_=rank_modp(J)
Jb=[J[k]+[(-b0[k])%p] for k in range(len(F))]
rJb,_=rank_modp(Jb)
print(f"FULL Jacobian mod p over {len(RF)} free vars:")
print(f"  rank(J)={rJ}, rank([J|-b])={rJb} -> {'CONSISTENT (Newton step exists mod p)' if rJ==rJb else 'INCONSISTENT mod p (no Newton step)'}")
# also how many eqs currently fail mod p
print("residuals mod p nonzero count:",sum(1 for x in b0 if x!=0),"/",len(F))
