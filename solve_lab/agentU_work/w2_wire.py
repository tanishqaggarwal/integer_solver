"""W2: what else constrains a leaf coordinate wire, and which coordinate carries m=1?"""
import pickle, collections
B='/home/user/integer_solver/solve_lab/agentU_work/'
D=pickle.load(open(B+'v_defs.pkl','rb')); A=pickle.load(open(B+'v_atoms.pkl','rb'))
AT=A['AT']; LEAFPIN=D['LEAFPIN']
L=pickle.load(open(B+'v_leaves.pkl','rb')); sel2exp=L['sel2exp']; p=L['p']; pts=L['pts']; shift=L['shift']
def vs(n,acc):
    if n[0]=='var': acc.add(n[1]); return
    if n[0]=='num': return
    for c in n[1:]: vs(c,acc)
def shape(n):
    k=n[0]
    if k=='var': return 'V'
    if k=='num': return 'C'
    return {'add':'(%s+%s)','sub':'(%s-%s)','mul':'(%s*%s)'}[k]%(shape(n[1]),shape(n[2]))
occ=collections.defaultdict(list)
for c,n in AT.items():
    s=set(); vs(n,s)
    for v in s: occ[v].append(c)
W={w for sel,w,C,z,m in LEAFPIN}
sh=collections.Counter()
for sel,w,C,z,m in LEAFPIN:
    for c in occ[w]:
        if AT[c][0]=='sub' and shape(AT[c]).startswith('((V*(V-C))'): continue   # its own pin
        sh[shape(AT[c])]+=1
print('shapes of NON-PIN atoms containing a leaf coordinate wire:')
for k,v in sh.most_common(): print('   %5d  %s'%(v,k))
print()
# which coordinate (X or Y in my model) carries m=1?
bysel=collections.defaultdict(dict)
for sel,w,C,z,m in LEAFPIN: bysel[sel][w]=(C%p,m,z)
tag={}
for s in sel2exp:
    X,Y=pts[s]; rawX=(X-shift)%p
    d={}
    for w,(C,m,z) in bysel[s].items():
        d['X' if C==rawX else ('Y' if C==Y else '?')]=(w,m,z)
    tag[s]=d
cnt=collections.Counter()
for s,d in tag.items():
    coord_m1=[k for k,(w,m,z) in d.items() if m==1]
    cnt[tuple(sorted(coord_m1))]+=1
print('which coordinate carries the m=1 pin, over the 256 leaves:', cnt.most_common())
mvals=collections.Counter(m for s,d in tag.items() for k,(w,m,z) in d.items() if m!=1)
print('distinct m>1 values:',len(mvals),'  sample:',list(mvals)[:6])
pickle.dump({'tag':tag}, open(B+'w_tag.pkl','wb'))
