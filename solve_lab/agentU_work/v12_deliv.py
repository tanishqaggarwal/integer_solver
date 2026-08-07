"""V12: end-to-end check of the decode against the 39,026 deliverable.
Reads the chord-gadget sites from Q (read-only) but all curve/exponent arithmetic is mine."""
import pickle, json, collections
B='/home/user/integer_solver/solve_lab/agentU_work/'
L=pickle.load(open(B+'v_leaves.pkl','rb'))
p=L['p']; shift=L['shift']; b=L['b']; sel2exp=L['sel2exp']; pts=L['pts']
Ndeg=115792089237316195423570985008687907852837564279074904382605163141518161494337
V=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
V={int(k[2:]) if k.startswith('x_') else int(k):int(v) for k,v in V.items()}
print('deliverable variables assigned:',len(V))
st=json.load(open('/home/user/integer_solver/solve_lab/agentQ_work/qstages.json'))['stages']
# wires hold raw u = X - shift ; add shift to land on Y^2 = X^3 + b
def pt(u,y):
    return ((V[u]+shift)%p, V[y]%p)
oncur=lambda P:(P[1]*P[1]-pow(P[0],3,p)-b)%p==0
# index every 2^e G by X-coordinate
e_of={}
for s,e in sel2exp.items(): e_of[pts[s][0]]=e
deg=[]; live=0; onc=0
for i,s in enumerate(st):
    try: A=pt(s['ua'],s['ya']); Bp=pt(s['ub'],s['yb'])
    except KeyError: continue
    if A==(shift%p,0) and Bp==(shift%p,0): continue
    live+=1
    if oncur(A) and oncur(Bp): onc+=1
    if A==Bp: deg.append((i,A))
print('stages with at least one non-zero input: %d/%d ; both inputs on the cubic: %d'%(live,len(st),onc))
print('stages whose two live inputs COINCIDE in the deliverable: %d'%len(deg))
for i,A in deg:
    print('   stage %d  x=%s...  is a leaf 2^e G ? e=%s'%(i,str(A[0])[:18],e_of.get(A[0])))
# the root: the stage with the largest live subtree, i.e. the one Q says has no parent
S=pickle.load(open(B+'v_supp2.pkl','rb')); supp=S['supp']; par=S['par']
def find(a):
    while par[a]!=a: a=par[a]
    return a
sizes=[(len(supp.get(find(s['u3']),())),i) for i,s in enumerate(st)]
sizes.sort(reverse=True)
print('largest u3 supports (size, stage):',sizes[:4])
