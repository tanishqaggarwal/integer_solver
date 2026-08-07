"""W1: every leaf pin is sel*(w-C) - m*z with z = a*b.  What are a and b?
If either factor is free, the leaf wire is free and the lie is FREE."""
import pickle, collections
B='/home/user/integer_solver/solve_lab/agentU_work/'
D=pickle.load(open(B+'v_defs.pkl','rb')); A=pickle.load(open(B+'v_atoms.pkl','rb'))
AT=A['AT']; EQ=A['EQ']; LEAFPIN=D['LEAFPIN']; DEFS=D['DEFS']; CONST=D['CONST']
L=pickle.load(open(B+'v_leaves.pkl','rb')); sel2exp=L['sel2exp']
def vs(n,acc):
    if n[0]=='var': acc.add(n[1]); return
    if n[0]=='num': return
    for c in n[1:]: vs(c,acc)
occ=collections.defaultdict(list)
for c,n in AT.items():
    s=set(); vs(n,s)
    for v in s: occ[v].append(c)
# z's defining product
prod={}
for v,lst in DEFS.items():
    for canon,sh,rv in lst:
        if sh=='(V-(V*V))':
            n=AT[canon]; prod[v]=(n[2][1][1], n[2][2][1])
free=lambda v: (v not in DEFS) and (v not in CONST)
rows=[]
for sel,w,C,z,m in LEAFPIN:
    ab=prod.get(z)
    if ab is None: rows.append((sel,w,z,m,None,None,None,None)); continue
    a,b=ab
    rows.append((sel,w,z,m,a,b,free(a),free(b)))
print('leaf pins:',len(rows))
print('z with a product definition:', sum(1 for r in rows if r[4] is not None))
print('pins where BOTH factors are free :', sum(1 for r in rows if r[6] and r[7]))
print('pins where AT LEAST ONE factor is free :', sum(1 for r in rows if r[6] or r[7]))
print('pins where NEITHER factor is free :', sum(1 for r in rows if r[4] is not None and not r[6] and not r[7]))
# classify the non-free factors
kinds=collections.Counter()
for sel,w,z,m,a,b,fa,fb in rows:
    if a is None: continue
    for v,f in ((a,fa),(b,fb)):
        if f: kinds['FREE']+=1
        elif v in CONST: kinds['const=%s'%CONST[v]]+=1
        else: kinds['defined:'+DEFS[v][0][1]]+=1
print('factor classification:',kinds.most_common())
# how many equations does each factor sit in?
fo=[]
for sel,w,z,m,a,b,fa,fb in rows:
    if a is None: continue
    for v,f in ((a,fa),(b,fb)):
        if f: fo.append((len(occ[v]),v,sel))
fo.sort()
print('free factors, atom-occurrence counts (lowest 8):', fo[:8])
pickle.dump({'rows':rows,'prod':prod}, open(B+'w_z.pkl','wb'))
