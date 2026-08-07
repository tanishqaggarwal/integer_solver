import pickle, collections
B='/home/user/integer_solver/solve_lab/agentU_work/'
T=pickle.load(open(B+'u_tree.pkl','rb')); D=pickle.load(open(B+'u_defs.pkl','rb'))
AT=pickle.load(open(B+'u_atoms.pkl','rb'))['ATOMS']
par=T['find']
def find(a):
    while par[a]!=a: a=par[a]
    return a
LIVE=T['LIVE']; supp=T['supp']
def vs(n,acc):
    if n[0]=='var': acc.add(n[1]); return
    if n[0]=='num': return
    for c in n[1:]: vs(c,acc)
byvar=collections.defaultdict(list)
for canon,n in AT.items():
    acc=set(); vs(n,acc)
    for v in acc: byvar[find(v)].append(canon)
for v in [22175]:
    print('=== node x%d support %d'%(v,len(supp[v])))
    for c in byvar[v][:40]: print('   ',c)
# how many live classes appear paired in some atom with another live class?
pairs=collections.Counter()
for canon,n in AT.items():
    acc=set(); vs(n,acc)
    lv=sorted({find(u) for u in acc} & set(LIVE))
    if len(lv)==2: pairs[tuple(lv)]+=1
print('atoms containing exactly 2 live classes:', len(pairs))
