"""Fast fold via LCA + cumulative coordinate swaps; enumerate small ON-sets against the target."""
import pickle, sys, time, itertools, collections
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
K=97553848499418123410591666447050222001188385549510401465815187079080512838891
M=pickle.load(open('/home/user/integer_solver/solve_lab/agentV_work/mirror/L/full_model.pkl','rb'))
C2=pickle.load(open('/home/user/integer_solver/solve_lab/agentV_work/mirror/L/calib2.pkl','rb'))
TGT=tuple(pickle.load(open('/home/user/integer_solver/solve_lab/agentV_work/mirror/L/target.pkl','rb')))
tree=M['tree']; PIN=M['PIN']; live=M['live']; ROOT=M['ROOT']; NODE=M['NODE']
perm=C2['perm']; ORIENT=C2['ORIENT']
parent={}; side_of={}
for n in NODE:
    for side,ch in (('va',NODE[n]['a']),('vb',NODE[n]['b'])): parent[ch]=n; side_of[ch]=side
swapup={x:(perm[(parent[x],side_of[x])][0]==1) for x in parent}
depth={ROOT:0}
def setdepth(n,d):
    depth[n]=d
    if tree[n] is not None:
        for c in tree[n]: setdepth(c,d+1)
setdepth(ROOT,0)
path={}; cums={}
for L in live:
    n=L; pl=[]; s=False; cl=[]
    while n!=ROOT:
        pl.append(n); cl.append(s); s^= swapup[n]; n=parent[n]
    pl.append(ROOT); cl.append(s)
    path[L]=pl; cums[L]=cl          # cums[L][i] = swap from L to pl[i]
ancset={L:set(path[L]) for L in live}
posn={L:{n:i for i,n in enumerate(path[L])} for L in live}
LEAF={L:tuple(PIN[L][1]) for L in live}
# swap from node n to ROOT
sw2root={}
for n in list(NODE)+list(tree):
    if n in sw2root: continue
    x=n; s=False
    while x!=ROOT: s^=swapup[x]; x=parent[x]
    sw2root[n]=s
def lca(A,B):
    sa=ancset[A]
    for n in path[B]:
        if n in sa: return n
    return ROOT
def sw(v,s): return (v[1],v[0]) if s else v
def chord(A,B,o):
    ax,ay,bx,by=A[o],A[1-o],B[o],B[1-o]
    d=(bx-ax)%p
    if d==0: return None
    l=(by-ay)*pow(d,p-2,p)%p
    ox=(l*l-ax-bx-K)%p; oy=(l*(ax-ox)-ay)%p
    return (ox,oy) if o==0 else (oy,ox)
def val_at(L,n):
    """leaf L's value expressed in node n's coord frame (n on L's path)"""
    return sw(LEAF[L],cums[L][posn[L][n]])
def fold2(A,B):
    m=lca(A,B)
    r=chord(val_at(A,m),val_at(B,m),ORIENT[m])
    return None if r is None else sw(r,sw2root[m])
if __name__=='__main__':
    t0=time.time()
    hits=[]
    # |S| = 1
    for L in live:
        if sw(LEAF[L],cums[L][-1])==TGT: hits.append((L,))
    print('|S|=1 done, hits',hits,' %.1fs'%(time.time()-t0))
    # |S| = 2
    t0=time.time(); n2=0
    for A,B in itertools.combinations(live,2):
        n2+=1
        if fold2(A,B)==TGT: hits.append((A,B)); print('HIT2',A,B)
    print('|S|=2 done (%d pairs) %.1fs'%(n2,time.time()-t0))
    pickle.dump(hits,open('hits.pkl','wb'))
