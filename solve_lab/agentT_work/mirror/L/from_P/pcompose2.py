#!/usr/bin/env python3
"""Agent P: composition test, hoisted.  cands[0] = smallest parent/child pair with two
live merges.  Distinguishes INDIVIDUAL lifting from SIMULTANEOUS lifting.
Every residual is recomputed directly from the shifted integers -- no expansion."""
import pickle,sys
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

# ---------- HOIST: every per-variable / per-block lookup done ONCE ----------
hvar={}          # var -> handle multiplier on its congruence atom (its legal lift step)
tri={}           # block -> [(c_k1, c_k2, c_k)]
for a in range(len(AP)):
    if pos[a]<19000: continue
    ap=AP[a]; vs={y for m in ap for y in m}
    hs=[x for x in vs if x in plift5.HANDLE]
    if not hs: continue
    ck=1
    for m,c in ap.items():
        if len(m)==1 and m[0] in plift5.HANDLE: ck=abs(c)
    for v in vs:
        if v not in plift5.HANDLE: hvar.setdefault(v,ck)
for j,b in enumerate(B):
    tri[j]=[(ca,cb,hvar.get(tc,1)) for (ca,cb,tc) in b['outs']]
print("hoist done: %d variable lift-steps, %d block condition triples"%(len(hvar),len(tri)))

def conds(j, iv):
    i1,i2,i3,i4,i5,i6=iv
    A=i1-i2; Bv=i4-i3; E=i1+i2+i5+Q
    N1=E*A*A-Bv*Bv; N2=A*(i3+i6)-Bv*(i2-i5)
    assert N1%P==0 and N2%P==0, "not mod-P valid"
    n1,n2=N1//P,N2//P
    return [(ca*n1+cb*n2, ck) for (ca,cb,ck) in tri[j]]

# ---------- cands[0] ----------
cands=[]
for jp,row in enumerate(F.SRC):
    for slot,k in enumerate(row):
        if k[0]!='S': continue
        jc=k[1]
        if not all(x[0]=='L' for x in F.SRC[jc]): continue
        if row[1-slot][0] not in ('L','S'): continue
        cands.append((len(F.supp[jp]),jp,jc,slot))
cands.sort()
print("candidates (support, parent, child, slot):",cands[:4])
_,jp,jc,slot=cands[0]
sel=set(F.supp[jp])
print("\n=== cands[0]: parent block %d <- child block %d , |S| = %d ==="%(jp,jc,len(sel)))
val,obs,und,sz,live=plift5.build(sel)
print("live blocks in this configuration:",sum(live.values()))

bc=[val[B[jc][k]] for k in ('i1','i2','i3','i4','i5','i6')]
bp=[val[B[jp][k]] for k in ('i1','i2','i3','i4','i5','i6')]
cc=conds(jc,bc); cp=conds(jp,bp)
print("child  conditions (residual mod c, c):",[(r%m,m) for r,m in cc])
print("parent conditions (residual mod c, c):",[(r%m,m) for r,m in cp])

# lift parameters
cpar=[hvar.get(B[jc][k],1) for k in ('i1','i2','i3','i4')]+[1,1]     # child: 4 leaves + i5,i6
cinp=(1,2) if slot==0 else (0,3)                                     # parent slots fed by child
ppar=[hvar.get(B[jp][k],1) for k in ('i1','i2','i3','i4')]+[1,1]
print("child  lift steps:",cpar)
print("parent lift steps:",ppar," (slots %s carry the child's output via a copy edge)"%(cinp,))

def solve(j, base, steps, idxs, label):
    """try to satisfy every non-vacuous condition of block j using parameters `idxs`."""
    allok=True
    for ridx,(r,m) in enumerate(conds(j,base)):
        if m==1: continue
        for q,e in sorted(fac(m).items()):
            qq=q**e; hit=None
            for pi in idxs:                       # single parameter sweep
                for t in range(qq):
                    iv=list(base); iv[pi]+=P*t*steps[pi]
                    if conds(j,iv)[ridx][0]%qq==0: hit=(pi,t); break
                if hit: break
            if hit is None:                       # small pair grid fallback
                for a_ in idxs:
                    for b_ in idxs:
                        if a_>=b_: continue
                        for t1 in range(min(qq,400)):
                            for t2 in range(min(qq,400)):
                                iv=list(base); iv[a_]+=P*t1*steps[a_]; iv[b_]+=P*t2*steps[b_]
                                if conds(j,iv)[ridx][0]%qq==0: hit=((a_,b_),(t1,t2)); break
                            if hit: break
                        if hit: break
                    if hit: break
            print("  %s cond %d  prime power %-9d -> %s"%(label,ridx,qq,hit))
            if hit is None: allok=False
    return allok

print("\n--- child conditions, using the child's own 6 parameters ---")
okc=solve(jc,bc,cpar,range(6),"child ")
print("--- parent conditions, using ONLY the copy-edge lifts (child untouched) ---")
okp=solve(jp,bp,ppar,list(cinp),"parent")
print("--- parent conditions, using all 6 parent parameters ---")
okp2=solve(jp,bp,ppar,range(6),"parent")
print()
print("=== VERDICT (cands[0], parent %d / child %d, |S|=%d) ==="%(jp,jc,len(sel)))
print("  child conditions solvable with child parameters      :",okc)
print("  parent conditions solvable with copy-edge lifts only :",okp)
print("  parent conditions solvable with all parent params    :",okp2)
print("  SIMULTANEOUS (child params and parent params disjoint):",okc and (okp or okp2))
