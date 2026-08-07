"""W3: for every slot, can two leaves in its two subtrees be driven to a COMMON
coordinate pair at zero pin cost?  Y is free (m=1 both sides).  X needs
  W = C_aX (mod M_a),  W = C_bX (mod M_b)   <=>  gcd(M_a,M_b) | (C_aX - C_bX)."""
import pickle, collections, math
B='/home/user/integer_solver/solve_lab/agentU_work/'
D=pickle.load(open(B+'v_defs.pkl','rb')); L=pickle.load(open(B+'v_leaves.pkl','rb'))
T=pickle.load(open(B+'v_tree_final.pkl','rb')); tg=pickle.load(open(B+'w_tag.pkl','rb'))['tag']
sel2exp=L['sel2exp']; exp2sel={v:k for k,v in sel2exp.items()}
rawC={}
for sel,w,C,z,m in D['LEAFPIN']: rawC[(sel,w)]=C
X={}   # exponent -> (wire, M, z, rawC)
Y={}
for s,d in tg.items():
    e=sel2exp[s]
    wx,mx,zx=d['X']; wy,my,zy=d['Y']
    X[e]=(wx,mx,zx,rawC[(s,wx)]); Y[e]=(wy,my,zy,rawC[(s,wy)])
print('m>1 always on X:',all(v[1]>1 for v in X.values()),'  m==1 always on Y:',all(v[1]==1 for v in Y.values()))
def feas(ea,eb):
    g=math.gcd(X[ea][1],X[eb][1])
    return (X[ea][3]-X[eb][3])%g==0, g
tot=0; ok=0; perslot={}
for k,(I,J) in T['children'].items():
    n=0; f=[]
    for ea in I:
        for eb in J:
            n+=1; good,g=feas(ea,eb)
            if good: f.append((ea,eb,g))
    tot+=n; ok+=len(f); perslot[frozenset((I,J))]=(len(I),len(J),n,len(f),f[:3])
print('slots: %d ; leaf pairs across a slot: %d ; ZERO-PIN-COST feasible pairs: %d'%(len(T['children']),tot,ok))
hit=[(v[0],v[1],v[3],v[4]) for v in perslot.values() if v[3]>0]
print('slots admitting at least one feasible pair: %d/%d'%(len(hit),len(T['children'])))
for h in sorted(hit,key=lambda z:-z[2])[:6]:
    print('   slot |I|=%d |J|=%d feasible=%d  e.g. %s'%h)
# distribution of gcds
gs=collections.Counter()
for ea in X:
    for eb in X:
        if ea<eb: gs[math.gcd(X[ea][1],X[eb][1])]+=1
print('gcd(M_a,M_b) distribution over all leaf pairs:',gs.most_common(8))
pickle.dump({'X':X,'Y':Y}, open(B+'w_xy.pkl','wb'))
