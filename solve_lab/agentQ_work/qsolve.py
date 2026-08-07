#!/usr/bin/env python3
"""Q-10b: general unit propagation over EVERY term of EVERY equation, mod p.

Each equation is a sum of terms, each of which must vanish.  Given the 256 selector bits and the
256 pinned leaf coordinates, repeatedly find a term with exactly one unknown variable and solve
for it (terms are linear or quadratic in any single variable).  No group theory enters: the chord
gadgets are solved as ordinary linear terms like everything else.  The root wire is then compared
against fold(S) computed independently from the ladder.
"""
import ast,re,json,sys,collections,pickle,os
p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
NV=38748
def const_val(n):
    if isinstance(n,ast.Constant): return n.value
    if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub) and isinstance(n.operand,ast.Constant): return -n.operand.value
    return None
def strip_outer(n):
    while True:
        if isinstance(n,ast.BinOp) and isinstance(n.op,ast.Mult):
            a,b=n.left,n.right; ca,cb=const_val(a),const_val(b)
            if ca is not None and cb is not None: return n
            if ca is not None: n=b; continue
            if cb is not None: n=a; continue
            if ast.unparse(a)==ast.unparse(b): n=a; continue
            return n
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub): n=n.operand; continue
        return n
def flatten(n):
    out=[]
    def rec(x):
        if isinstance(x,ast.BinOp) and isinstance(x.op,ast.Add): rec(x.left); rec(x.right)
        else: out.append(x)
    rec(n); return out
CACHE='qterms.pkl'
if os.path.exists(CACHE):
    terms=pickle.load(open(CACHE,'rb'))
else:
    VAR=re.compile(r'x_(\d+)')
    seen={}; terms=[]
    for L in open('/home/user/integer_solver/EQUATIONS.txt'):
        L=L.strip()
        if not L: continue
        node=ast.parse(L.rsplit('=',1)[0].strip(),mode='eval').body
        for t in flatten(strip_outer(node)):
            t=strip_outer(t)
            s=ast.unparse(t)
            if s in seen: continue
            vs=tuple(sorted({int(m) for m in VAR.findall(s)}))
            if not vs: continue
            seen[s]=1; terms.append((s,vs))
    pickle.dump(terms,open(CACHE,'wb'))
print('distinct terms:',len(terms))
byvar=collections.defaultdict(list)
for i,(s,vs) in enumerate(terms):
    for v in vs: byvar[v].append(i)
code=[None]*len(terms)
VARSUB=re.compile(r'x_(\d+)')
def getcode(i):
    if code[i] is None:
        code[i]=compile(VARSUB.sub(r'V[\1]',terms[i][0]),'<t>','eval')
    return code[i]
def propagate(V):
    ns={'V':V,'__builtins__':{}}
    q=collections.deque(v for v in range(NV) if V[v] is not None)
    inq=set(q); solved=0; contra=0
    while q:
        v=q.popleft(); inq.discard(v)
        for i in byvar[v]:
            unk=[u for u in terms[i][1] if V[u] is None]
            if len(unk)!=1: continue
            u=unk[0]; c=getcode(i)
            f=[]
            for tv in (0,1,2):
                V[u]=tv; f.append(eval(c,ns)%p)
            V[u]=None
            if (f[2]-2*f[1]+f[0])%p==0:
                a=(f[1]-f[0])%p; b=f[0]%p
                if a==0:
                    if b!=0: contra+=1
                    continue
                V[u]=(-b)*pow(a,p-2,p)%p
            else:
                A=((f[2]-2*f[1]+f[0])*pow(2,p-2,p))%p
                B=((4*f[1]-3*f[0]-f[2])*pow(2,p-2,p))%p
                C=f[0]%p
                D=(B*B-4*A*C)%p
                r=pow(D,(p+1)//4,p)
                if r*r%p!=D: continue
                if D==0: V[u]=(-B)*pow(2*A,p-2,p)%p
                else: continue
            solved+=1
            if u not in inq: q.append(u); inq.add(u)
    return solved,contra
if __name__=='__main__':
    import random
    sys.path.insert(0,'.')
    from qgrp import add,p as _p,cs
    ROOTX,ROOTY=24468,18956
    leaf={int(g):v for g,v in json.load(open('qleaf.json')).items()}
    LP={g:(int(v[0]),int(v[1])) for g,v in leaf.items()}
    LW={g:(v[2],v[3]) for g,v in leaf.items()}
    lad=json.load(open('qladder.json')); e2s={int(k):v for k,v in lad['exp2sel'].items()}
    def fold(S):
        F=None
        for g in S: F=add(F,LP[g])
        return F
    random.seed(3); sel=sorted(LP)
    sizes=[int(x) for x in (sys.argv[1:] or ['1','2','5','17','64','128','200','256'])]
    print('%-5s %-9s %-9s %-9s %s'%('|S|','solved','rootX','rootY','fold(S) match'))
    for n in sizes:
        S=set(random.sample(sel,n)) if n<256 else set(sel)
        V=[None]*NV
        for g in LP:
            V[g]=1 if g in S else 0
            wx,wy=LW[g]; V[wx]=(LP[g][0]-cs)%p; V[wy]=LP[g][1]%p
        s,c=propagate(V)
        F=fold(S)
        wantx=(F[0]-cs)%p if F else None; wanty=F[1]%p if F else None
        m = (V[ROOTX]==wantx and V[ROOTY]==wanty) if F else None
        print('%-5d %-9d %-9s %-9s %s'%(n,s,V[ROOTX] is not None,V[ROOTY] is not None,m))
