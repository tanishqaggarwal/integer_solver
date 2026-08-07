import sys, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentO_work')
import simO, engine as E, harness as H, full11 as F, fast, sparse
d=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
vd=[0]*E.NV
for k,x in d.items(): vd[int(k.split('_')[1])]=int(x)
r=F.solve_pair(24601,2081,verbose=False); n,ns,av=r
v0=E.forward(ns); bad0=E.badatoms(v0)
print('pair state bad0',sorted(bad0))
CARR=[8731,9118,9413,1329,10903,17325]
ROWS=[20649,20652,32148]
CF=set()
for a in ROWS: CF|=set(E.cone(a)[1])
print('cone free vars of trio:',len(CF),'carriers in cone:',[c for c in CARR if c in CF])
for f in CARR:
    b1,_=fast.resid_delta(v0,bad0,{f:v0[f]+1})
    b2,_=fast.resid_delta(v0,bad0,{f:v0[f]+2})
    col={a:b1.get(a,0)-bad0.get(a,0) for a in set(b1)|set(bad0) if b1.get(a,0)-bad0.get(a,0)}
    nl=[a for a in set(b2)|set(bad0)|set(col) if b2.get(a,0)-bad0.get(a,0)!=2*col.get(a,0)]
    print(f'x_{f}: touches {len(col)} atoms, on trio: '
          f'{[(a,len(str(abs(col.get(a,0))))) for a in ROWS]} nonlin_on_trio={[a for a in ROWS if a in nl]} tot_nonlin={len(nl)}')
# now try the deliverable's carrier values
s2=dict(ns)
for f in CARR: s2[f]=vd[f]
v2=E.forward(s2); a2=E.badatoms(v2)
print('with deliverable carrier values: bad=',sorted(a2),'fails=',len(E.eqfails(a2)))
