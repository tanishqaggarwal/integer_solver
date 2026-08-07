#!/usr/bin/env python3
"""Agent P: composition test.  Two live merges in a parent/child relation.
NO symbolic expansion anywhere -- every residual is recomputed directly from the
shifted integers, which is the guard that caught the A*h2/B*h2 bug last round."""
import pickle,sys
from math import gcd
sys.set_int_max_str_digits(10**7)
W='/home/user/integer_solver/solve_lab/agentP_work/'
D=pickle.load(open(W+'model4.pkl','rb')); AP=D['AP']
S=pickle.load(open(W+'slp.pkl','rb')); topo=S['topo']; outof=S['outof']
B=pickle.load(open(W+'blocks.pkl','rb'))
import plift5, pfold as F
P=plift5.P; Q=F.Q
pos={a:i for i,a in enumerate(topo)}

def fac(n):
    f={};d=2
    while d*d<=n:
        while n%d==0: f[d]=f.get(d,0)+1; n//=d
        d+=1
    if n>1: f[n]=f.get(n,0)+1
    return f

def stepof(v):
    """legal lift step for variable v, as a multiple of P.  Returns (step, found)."""
    for a in range(len(AP)):
        if pos[a]<19000: continue
        ap=AP[a]; vs={y for m in ap for y in m}
        if v in vs and any(x in plift5.HANDLE for x in vs):
            for m,c in ap.items():
                if len(m)==1 and m[0] in plift5.HANDLE: return abs(c),True
    return 1,False

# ---- 1. verify the i3/i4 step lookup (flagged last round) ----
for j,row in enumerate(F.SRC):
    if all(k[0]=='L' for k in row): jj=j; break
b2=B[jj]
print("i3/i4 step audit at block %d:"%jj)
for nm in ('i1','i2','i3','i4'):
    s,found=stepof(b2[nm])
    print("   %s = x%-6d  step=%-9d  atom found=%s"%(nm,b2[nm],s,found))

# ---- 2. find a parent/child pair of live merges ----
pair=None
for jp,row in enumerate(F.SRC):
    kinds=[k[0] for k in row]
    if 'S' not in kinds: continue
    for slot,k in enumerate(row):
        if k[0]!='S': continue
        jc=k[1]
        if all(x[0]=='L' for x in F.SRC[jc]) and row[1-slot][0]=='L':
            pair=(jp,jc,slot); break
    if pair: break
jp,jc,slot=pair
sel=set()
for x in F.SRC[jc]: sel.add(x[1])
sel.add(F.SRC[jp][1-slot][1])
print("\nparent block %d  <-  child block %d ; |S| = %d"%(jp,jc,len(sel)))

val,obs,und,sz,live=plift5.build(sel)
print("live blocks in this configuration:",sum(live.values()))

# ---- 3. the copy edge that carries the child's output into the parent ----
# parent input pair on `slot` should be the child's output, joined by a mod-P copy
cinp=(B[jp]['i2'],B[jp]['i3']) if slot==0 else (B[jp]['i1'],B[jp]['i4'])
print("parent inputs on the child slot: x%d , x%d"%cinp)
for v in cinp:
    s,found=stepof(v)
    print("   x%-6d copy-edge lift step = %-9d (atom found=%s)"%(v,s,found))

# ---- 4. residuals, by DIRECT recomputation ----
def conds(blk, iv):
    """iv = the six integer inputs; return list of (residual_over_P, modulus)."""
    i1,i2,i3,i4,i5,i6=iv
    A=i1-i2; Bv=i4-i3; E=i1+i2+i5+Q
    N1=E*A*A-Bv*Bv; N2=A*(i3+i6)-Bv*(i2-i5)
    assert N1%P==0 and N2%P==0, "base not mod-P valid"
    n1,n2=N1//P,N2//P
    out=[]
    for (ca,cb,tc) in blk['outs']:
        cg=[a for a in range(len(AP)) if tc in {y for m in AP[a] for y in m} and pos[a]>19000]
        ck=1
        for m,c in AP[cg[0]].items():
            if len(m)==1 and m[0] in plift5.HANDLE: ck=abs(c)
        out.append((ca*n1+cb*n2, ck))
    return out

base_c=[val[B[jc][k]] for k in ('i1','i2','i3','i4','i5','i6')]
base_p=[val[B[jp][k]] for k in ('i1','i2','i3','i4','i5','i6')]
cc=conds(B[jc],base_c); cp=conds(B[jp],base_p)
print("\nchild  conditions (residual mod c, modulus c):",[(r%m,m) for r,m in cc])
print("parent conditions (residual mod c, modulus c):",[(r%m,m) for r,m in cp])
nontriv=[m for _,m in cc+cp if m>1]
print("non-vacuous conditions in this configuration:",len(nontriv),nontriv)

# ---- 5. do the parent's conditions decouple via the copy-edge lifts? ----
print("\n=== COMPOSITION TEST ===")
csteps=[stepof(v)[0] for v in cinp]
print("copy-edge lift steps available to the parent:",csteps)
ok=True
for ridx,(r,m) in enumerate(cp):
    if m==1: print("  parent cond %d: modulus 1, vacuous"%ridx); continue
    solved={}
    for q,e in sorted(fac(m).items()):
        qq=q**e; hit=None
        for u1 in range(qq):
            for u2 in range(qq):
                iv=list(base_p)
                if slot==0: iv[1]+=P*u1*csteps[0]; iv[2]+=P*u2*csteps[1]
                else:       iv[0]+=P*u1*csteps[0]; iv[3]+=P*u2*csteps[1]
                rr=conds(B[jp],iv)[ridx][0]
                if rr%qq==0: hit=(u1,u2); break
            if hit: break
        solved[qq]=hit
        print("  parent cond %d, prime power %-8d -> %s"%(ridx,qq,hit))
        if hit is None: ok=False
print("\nparent conditions solvable using ONLY the copy-edge lifts (child untouched):",ok)
