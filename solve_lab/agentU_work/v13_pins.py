"""V13: what actually pins a leaf coordinate wire?  The pin is  sel*(w - C) - m*z = 0,
so with sel=1 the wire carries C + m*z.  Measure z's other constraints, in EQUATIONS."""
import pickle, collections
B='/home/user/integer_solver/solve_lab/agentU_work/'
D=pickle.load(open(B+'v_defs.pkl','rb')); A=pickle.load(open(B+'v_atoms.pkl','rb'))
AT=A['AT']; EQ=A['EQ']; LEAFPIN=D['LEAFPIN']; CONST=D['CONST']; DEFS=D['DEFS']
L=pickle.load(open(B+'v_leaves.pkl','rb')); sel2exp=L['sel2exp']
def vs(n,acc):
    if n[0]=='var': acc.add(n[1]); return
    if n[0]=='num': return
    for c in n[1:]: vs(c,acc)
atomvars={c:(lambda a:(vs(a,s) or s))(n) for c,n in ((c,n) for c,n in AT.items()) for s in [set()]}
occ=collections.defaultdict(list)
for c,vv in atomvars.items():
    for v in vv: occ[v].append(c)
# equation footprint per atom
eqof=collections.defaultdict(set)
for i,lst in enumerate(EQ):
    for c in lst: eqof[c].add(i)
print('leaf pins: %d ; multipliers: %s'%(len(LEAFPIN), collections.Counter('m=1' if m==1 else 'm>1' for *_,m in LEAFPIN)))
rows=[]
for sel,w,C,z,m in LEAFPIN:
    zc=occ[z]; wc=occ[w]
    zdef = z in DEFS or z in CONST
    rows.append((sel,w,z,m,len(zc),len(wc),zdef, z in CONST, CONST.get(z)))
print('z wires that occur in exactly 1 atom (the pin itself):', sum(1 for r in rows if r[4]==1))
print('z wires with a definition or constant pin:', sum(1 for r in rows if r[6]))
print('z occurrence-count histogram:', collections.Counter(r[4] for r in rows).most_common())
print('z pinned to a constant:', collections.Counter(str(r[8]) for r in rows if r[7]).most_common(5))
print('w occurrence-count histogram:', collections.Counter(r[5] for r in rows).most_common(8))
# what else contains z?
ex=collections.Counter()
def shape(n):
    k=n[0]
    if k=='var': return 'V'
    if k=='num': return 'C'
    return {'add':'(%s+%s)','sub':'(%s-%s)','mul':'(%s*%s)'}[k]%(shape(n[1]),shape(n[2]))
for sel,w,C,z,m in LEAFPIN:
    for c in occ[z]:
        ex[shape(AT[c])]+=1
print('shapes of atoms containing a z wire:', ex.most_common(10))
# --- the price: equations touched by the pin atom of each leaf, and by z ---
pinatom={}
for c,n in AT.items():
    sh=shape(n)
    if sh=='((V*(V-C))-V)': pinatom[(n[1][1][1],n[1][2][1][1])]=c
    elif sh=='((V*(V-C))-(C*V))': pinatom[(n[1][1][1],n[1][2][1][1])]=c
costs=[]
for sel,w,C,z,m in LEAFPIN:
    c=pinatom[(sel,w)]
    E=set(eqof[c])
    for c2 in occ[z]: E|=eqof[c2]
    costs.append((len(E),sel,w,z,m))
costs.sort()
print('equations touched by {pin atom} U {every atom containing z}: min %d  median %d  max %d'%(
      costs[0][0], costs[len(costs)//2][0], costs[-1][0]))
print('cheapest 6:',[(c[0],'sel x%d exp %d'%(c[1],sel2exp[c[1]]),'m=%d'%c[4]) for c in costs[:6]])
pickle.dump({'rows':rows,'costs':costs}, open(B+'v_pins.pkl','wb'))
