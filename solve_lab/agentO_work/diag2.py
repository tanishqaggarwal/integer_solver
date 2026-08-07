import sys, time, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentO_work')
import simO, engine as E, sparse
bit=int(sys.argv[1]) if len(sys.argv)>1 else 22492
s=dict(simO.C.base)
if bit: s[bit]=1
v0=E.forward(s); bad0=E.badatoms(v0)
S,cols,nonlin,rounds=simO.closure(v0,bad0,{bit},5,4000)
print('knobs',len(S),'bad0',sorted(bad0))
nlf={}
for f,a in nonlin: nlf.setdefault(f,set()).add(a)
clean=[f for f in S if f not in nlf]
print('fully-linear knobs',len(clean))
# per-row: how many knobs reach it, how many linearly
for a in sorted(bad0):
    reach=[f for f in S if a in cols[f]]
    lin=[f for f in reach if (f,a) not in nonlin]
    cl=[f for f in lin if f in set(clean)]
    print(f'  row {a}: reach={len(reach)} lin={len(lin)} cleanknob={len(cl)} rhs_bits={abs(bad0[a]).bit_length()}')
