import pickle, collections, sys
B='/home/user/integer_solver/solve_lab/agentU_work/'
D=pickle.load(open(B+'u_atoms.pkl','rb'))
A=D['ATOMS']
def shape(n):
    k=n[0]
    if k=='var': return 'V'
    if k=='num': return 'C'
    if k=='neg': return '-'+shape(n[1])
    return {'add':'(%s+%s)','sub':'(%s-%s)','mul':'(%s*%s)'}[k]%(shape(n[1]),shape(n[2]))
def vs(n,acc):
    if n[0]=='var': acc.add(n[1]); return
    if n[0]=='num': return
    for c in n[1:]: vs(c,acc)

DEFS=collections.defaultdict(list)   # var -> list of (canon, shape, rhs_vars)
COPY=[]                              # (a,b) pure copy
CONST={}                             # var -> constant
LEAFPIN=[]                           # (sel, w, C, z)
OTHER=collections.Counter()
allvars=set()
for canon,n in A.items():
    sh=shape(n); acc=set(); vs(n,acc); allvars|=acc
    if sh=='(V-V)': COPY.append((n[1][1],n[2][1]))
    elif sh=='(V-C)': CONST.setdefault(n[1][1],n[2][1])
    elif sh=='(C-V)': CONST.setdefault(n[2][1],n[1][1])
    elif sh.startswith('(V-') and n[1][0]=='var':
        rv=set(); vs(n[2],rv); DEFS[n[1][1]].append((canon,sh,rv))
    elif sh=='((V*(V-C))-V)':
        # (sel*(w-C)) - z
        LEAFPIN.append((n[1][1][1], n[1][2][1][1], n[1][2][2][1], n[2][1], 1))
    elif sh=='((V*(V-C))-(C*V))':
        LEAFPIN.append((n[1][1][1], n[1][2][1][1], n[1][2][2][1], n[2][2][1], n[2][1][1]))
    else: OTHER[sh]+=1
print('vars seen', len(allvars), 'min',min(allvars),'max',max(allvars))
print('DEFS vars', len(DEFS), ' multi-def vars', sum(1 for v in DEFS.values() if len(v)>1))
print('COPY', len(COPY), 'CONST', len(CONST), 'LEAFPIN', len(LEAFPIN))
print('OTHER shapes', OTHER.most_common())
sels=collections.Counter(l[0] for l in LEAFPIN)
print('distinct sel vars in leaf pins:', len(sels), 'counts', collections.Counter(sels.values()))
ws=collections.Counter(l[1] for l in LEAFPIN); print('distinct w wires', len(ws))
Cs=set(l[2] for l in LEAFPIN); print('distinct big constants', len(Cs), 'bitlens', sorted(set(c.bit_length() for c in Cs))[:5],'..',sorted(set(c.bit_length() for c in Cs))[-3:])
pickle.dump({'DEFS':dict(DEFS),'COPY':COPY,'CONST':CONST,'LEAFPIN':LEAFPIN,'allvars':allvars}, open(B+'u_defs.pkl','wb'))
