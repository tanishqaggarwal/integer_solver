"""Standalone fast fold evaluator over the calibrated tree (no Engine needed)."""
import pickle, sys, collections
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
K=97553848499418123410591666447050222001188385549510401465815187079080512838891
M=pickle.load(open('/home/user/integer_solver/solve_lab/agentT_work/mirror/L/full_model.pkl','rb'))
C2=pickle.load(open('/home/user/integer_solver/solve_lab/agentT_work/mirror/L/calib2.pkl','rb'))
tree=M['tree']; PIN=M['PIN']; live=M['live']; ROOT=M['ROOT']; order=M['order']; sub=M['sub']
perm=C2['perm']; ORIENT=C2['ORIENT']; NODE=M['NODE']
parent={}
for n in NODE:
    for side,ch in (('va',NODE[n]['a']),('vb',NODE[n]['b'])): parent[ch]=(n,side)
def chord(A,B,o):
    ax,ay,bx,by=A[o],A[1-o],B[o],B[1-o]
    d=(bx-ax)%p
    if d==0: return None
    l=(by-ay)*pow(d,p-2,p)%p
    ox=(l*l-ax-bx-K)%p; oy=(l*(ax-ox)-ay)%p
    return (ox,oy) if o==0 else (oy,ox)
def invchord(O,A,o):
    """B with chord(A,B,o)==O"""
    ax,ay,ox,oy=A[o],A[1-o],O[o],O[1-o]
    d=(ax-ox)%p
    if d==0: return None
    l=(oy+ay)*pow(d,p-2,p)%p
    bx=(l*l-ax-ox-K)%p; by=(ay+l*(bx-ax))%p
    return (bx,by) if o==0 else (by,bx)
LEAFVAL={L:tuple(PIN[L][1]) for L in live}
def fold(S):
    S=set(S); val={}; isl={}
    for L in tree:
        if tree[L] is None: isl[L]=L in S; val[L]=LEAFVAL.get(L) if L in S else None
    for n in order:
        a,b=tree[n]; la,lb=isl[a],isl[b]; isl[n]=la or lb
        def pr(ch,side):
            pm=perm[(n,side)]; v=val[ch]
            return None if v is None else (v[pm[0]],v[pm[1]])
        if la and lb:
            o=ORIENT[n]; A=pr(a,'va'); B=pr(b,'vb')
            val[n]=None if (A is None or B is None or o=='DEAD') else chord(A,B,o)
        elif la: val[n]=pr(a,'va')
        elif lb: val[n]=pr(b,'vb')
        else: val[n]=None
    return val[ROOT]
# sparse fold: only nodes on paths from S to root
anc={}
def ancestors(L):
    if L in anc: return anc[L]
    r=[]; n=L
    while n!=ROOT:
        n,side=parent[n]; r.append(n)
    anc[L]=r; return r
def fold_sparse(S):
    S=set(S)
    nodes=set()
    for L in S: nodes.update(ancestors(L))
    nl=[n for n in order if n in nodes]
    val={}; isl={}
    for L in S: val[L]=LEAFVAL[L]
    for n in nl:
        a,b=tree[n]
        A=val.get(a); B=val.get(b)
        def pr(v,side):
            pm=perm[(n,side)]; return None if v is None else (v[pm[0]],v[pm[1]])
        A=pr(A,'va'); B=pr(B,'vb')
        if A is not None and B is not None:
            o=ORIENT[n]; val[n]=chord(A,B,o) if o!='DEAD' else None
        else: val[n]=A if A is not None else B
    return val.get(ROOT)
if __name__=='__main__':
    import random,time
    rnd=random.Random(5)
    for k in (1,2,5,40):
        S=rnd.sample(live,k)
        a=fold(S); b=fold_sparse(S)
        print(k,'match',a==b, a)
    t0=time.time()
    for _ in range(2000): fold_sparse(rnd.sample(live,3))
    print('fold_sparse |S|=3: %.1f us'%((time.time()-t0)/2000*1e6))
