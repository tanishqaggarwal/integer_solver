import sys, time, json, math
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentO_work')
import simO, engine as E, sparse, harness as H
bit=int(sys.argv[1]) if len(sys.argv)>1 else 22492
s=dict(simO.C.base)
if bit: s[bit]=1
v0=E.forward(s); bad0=E.badatoms(v0)
S,cols,nonlin,rounds=simO.closure(v0,bad0,{bit} if bit else set(),6,8000)
print('knobs',len(S))
P=115792089237316195423570985008687907853269984665640564039457584007908834671663
for a in sorted(bad0):
    reach={f:cols[f][a] for f in S if a in cols[f]}
    lin={f:c for f,c in reach.items() if (f,a) not in nonlin}
    g=0
    for c in lin.values(): g=math.gcd(g,c)
    r=bad0[a]
    print(f'row {a}: atom="{H.atoms[a][:90]}"')
    print(f'   reach={len(reach)} lin={len(lin)} gcd(lin)={g if g<10**12 else str(g)[:20]+"..("+str(len(str(g)))+"d)"} '
          f'rhs%gcd={"0" if g and r%g==0 else "NZ"} solo_solvable={bool(g) and r%g==0}')
    for f,c in sorted(lin.items())[:8]:
        print(f'      knob x_{f}: coef {str(c)[:40]}..({len(str(abs(c)))}d)  p|c={c%P==0}')
