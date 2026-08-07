import sys, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentO_work')
import simO, engine as E, harness as H
d=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
vd=[0]*E.NV
for k,x in d.items(): vd[int(k.split('_')[1])]=int(x)
FR=[u for u in range(E.NV) if E.definer[u] is None]
seed={u:vd[u] for u in FR if vd[u]!=0}
v=E.forward(seed)
pos={u:k for k,u in enumerate(E.SEQ)}
mism=sorted([u for u in range(E.NV) if v[u]!=vd[u]], key=lambda u: pos.get(u,-1))
print('mismatches',len(mism))
for u in mism:
    i,kind=E.definer[u]
    ins=[w for w in H.avars[i] if w!=u]
    bad_in=[w for w in ins if v[w]!=vd[w]]
    print(f'seq{pos[u]:>6} x_{u} kind={kind} atom={H.atoms[i][:80]}')
    print(f'        inputs_mismatched={bad_in}  ROOT={not bad_in}')
