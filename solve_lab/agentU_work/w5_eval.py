"""W5: forward-only evaluator built from MY parse.  Every defined variable is recomputed
from its OWN definition, in SLP order.  Nothing is ever solved backwards.
Calibration: propagating the deliverable must reproduce the deliverable's score."""
import sys, json, pickle, collections, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab')
B='/home/user/integer_solver/solve_lab/agentU_work/'
D=pickle.load(open(B+'v_defs.pkl','rb')); A=pickle.load(open(B+'v_atoms.pkl','rb'))['AT']
DEFS=D['DEFS']; COPY=D['COPY']; CONST=D['CONST']; NV=38748
OPS={}   # var -> (kind, args)
for v,lst in DEFS.items():
    canon,sh,rv=lst[0]; n=A[canon]; r=n[2]
    if   sh=='(V-(V*V))': OPS[v]=('mul',r[1][1],r[2][1])
    elif sh=='(V-(V+V))': OPS[v]=('add',r[1][1],r[2][1])
    elif sh=='(V-(V-V))': OPS[v]=('sub',r[1][1],r[2][1])
    elif sh=='(V-(C*V))': OPS[v]=('cmul',r[1][1],r[2][1])
    elif sh=='(V-(V*C))': OPS[v]=('cmul',r[2][1],r[1][1])
    elif sh=='(V-(V+C))': OPS[v]=('cadd',r[2][1],r[1][1])
    elif sh=='(V-(C-V))': OPS[v]=('csub',r[1][1],r[2][1])
for v,c in CONST.items(): OPS.setdefault(v,('const',c,None))
# copies: orient towards the side that already has a definition
cp=[]
for a,b in COPY:
    if a in OPS and b not in OPS: cp.append((b,a))
    elif b in OPS and a not in OPS: cp.append((a,b))
    else: cp.append((a,b))
for t,s in cp: OPS.setdefault(t,('copy',s,None))
print('defined variables:',len(OPS),' free:',NV-len(OPS))
# topological order (Kahn); anything left over is in a cycle and is evaluated last, iteratively
deps={v:[x for x in (o[1],o[2]) if isinstance(x,int) and o[0]!='const' and not (o[0] in('cmul','cadd','csub') and x==o[1] and False)] for v,o in OPS.items()}
for v,o in OPS.items():
    if o[0]=='const': deps[v]=[]
    elif o[0]=='copy': deps[v]=[o[1]]
    elif o[0] in ('cmul','cadd'): deps[v]=[o[1]]
    elif o[0]=='csub': deps[v]=[o[2]]
    else: deps[v]=[o[1],o[2]]
indeg=collections.Counter(); users=collections.defaultdict(list)
for v,ds in deps.items():
    for u in ds:
        if u in OPS: indeg[v]+=1; users[u].append(v)
q=collections.deque([v for v in OPS if indeg[v]==0]); order=[]
while q:
    v=q.popleft(); order.append(v)
    for w in users.get(v,()):
        indeg[w]-=1
        if indeg[w]==0: q.append(w)
cyc=[v for v in OPS if v not in set(order)]
print('topological order length %d ; variables in cycles: %d'%(len(order),len(cyc)))
ORDER=order+cyc
def propagate(v, rounds=3):
    v=list(v)
    for _ in range(rounds):
        for t in ORDER:
            k,a,b=OPS[t]
            if   k=='const': v[t]=a
            elif k=='copy':  v[t]=v[a]
            elif k=='mul':   v[t]=v[a]*v[b]
            elif k=='add':   v[t]=v[a]+v[b]
            elif k=='sub':   v[t]=v[a]-v[b]
            elif k=='cmul':  v[t]=a*v[b]
            elif k=='cadd':  v[t]=a+v[b]
            elif k=='csub':  v[t]=a-v[b]
        if not cyc: break
    return v
pickle.dump({'OPS':OPS,'ORDER':ORDER,'cyc':cyc}, open(B+'w_eval.pkl','wb'))
if __name__=='__main__':
    import checker
    codes,_=checker.load_equations()
    v0=checker.load_assignment('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
    print('deliverable          : %d failing'%len(checker.evaluate_all(codes,v0)))
    v1=propagate(v0)
    diff=sum(1 for i in range(NV) if v0[i]!=v1[i])
    print('propagated deliverable: %d failing ; variables changed by propagation: %d'%(
        len(checker.evaluate_all(codes,v1)), diff))
