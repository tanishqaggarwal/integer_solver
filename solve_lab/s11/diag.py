import sys, os, json, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, engine
from zsolve import solve_int
P=L.P
CTRL=json.load(open('controls.json'))
BAD=[26719,26721,26723,26733,28438,32342,36185]
def resid(theta):
    v=engine.apply_theta(theta)
    return [fw.evalpoly(L.polys[a],v) for a in BAD]
theta={c:0 for c in CTRL}
r0=resid(theta)
J=[[0]*len(CTRL) for _ in BAD]
for j,c in enumerate(CTRL):
    th=dict(theta); th[c]=1
    r1=resid(th)
    for i in range(len(BAD)): J[i][j]=r1[i]-r0[i]
json.dump({'r0':[str(x) for x in r0],'J':[[str(x) for x in row] for row in J],'CTRL':CTRL}, open('jac.json','w'))
print("private handle deltas mod P:")
for i,a in enumerate(BAD):
    for j,c in enumerate(CTRL):
        if J[i][j] and len(L.var_atoms[c])==1:
            d=J[i][j]
            print(f"  a{a} x{c}: delta/P = {d//P if d%P==0 else 'NOT mult of P'}   delta%P={d%P!=0}")
print()
print("row gcd vs residual (integer solvability per row):")
for i,a in enumerate(BAD):
    g=0
    for j in range(len(CTRL)): g=math.gcd(g,J[i][j])
    print(f"  a{a}: gcd={str(g)[:30]}({len(str(g))}d)  r0%gcd==0: {r0[i]%g==0}")
print()
print("mod-P system (rows: mirror checks, unknowns x5096,x19750):")
i5=CTRL.index(5096); i19=CTRL.index(19750)
rows=[]
for i,a in enumerate(BAD):
    rows.append((a, r0[i]%P, J[i][i5]%P, J[i][i19]%P))
    print(f"  a{a}: {J[i][i5]%P!=0} {J[i][i19]%P!=0}  r0modP={'0' if r0[i]%P==0 else 'nz'}")
