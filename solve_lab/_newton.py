import heal_harness as H, pickle, time, random
p=H.p
AD=pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/adcache.pkl','rb'))
gate_parsed=AD['gate_parsed']
def inv(x): return pow(x%p,p-2,p)
class Dm:
    __slots__=('a','b')
    def __init__(s,a,b=0): s.a=a%p; s.b=b%p
    def __add__(s,o):
        o=o if isinstance(o,Dm) else Dm(o); return Dm(s.a+o.a,s.b+o.b)
    __radd__=__add__
    def __sub__(s,o):
        o=o if isinstance(o,Dm) else Dm(o); return Dm(s.a-o.a,s.b-o.b)
    def __rsub__(s,o):
        o=o if isinstance(o,Dm) else Dm(o); return Dm(o.a-s.a,o.b-s.b)
    def __mul__(s,o):
        o=o if isinstance(o,Dm) else Dm(o); return Dm(s.a*o.a,s.a*o.b+s.b*o.a)
    __rmul__=__mul__
    def __neg__(s): return Dm(-s.a,-s.b)
def tangent(seedmap,val):
    dv=[0]*H.NVARS
    for f,c in seedmap.items(): dv[f]=c%p
    for t in H.order:
        pr=gate_parsed[t]; op=pr[0]
        if op=='id':
            k=pr[1]; dv[t]=dv[k[1]] if k[0]=='v' else 0
        elif op=='add':
            k1,k2=pr[1],pr[2]; dv[t]=((dv[k1[1]] if k1[0]=='v' else 0)+(dv[k2[1]] if k2[0]=='v' else 0))%p
        elif op=='sub':
            k1,k2=pr[1],pr[2]; dv[t]=((dv[k1[1]] if k1[0]=='v' else 0)-(dv[k2[1]] if k2[0]=='v' else 0))%p
        elif op=='mul':
            k1,k2=pr[1],pr[2]
            v1=val[k1[1]] if k1[0]=='v' else k1[1]; v2=val[k2[1]] if k2[0]=='v' else k2[1]
            d1=dv[k1[1]] if k1[0]=='v' else 0; d2=dv[k2[1]] if k2[0]=='v' else 0
            dv[t]=(v1*d2+v2*d1)%p
    return dv
def build_and_solve(val, Fset):
    # closure iter-2 for robustness
    gS=tangent  # not used
    # movers of S,T
    import copy
    # gradient movers: reuse: frees with nonzero tangent effect on S or T -> approximate by initial R via anc
    R=set(H.anc[35389]|H.anc[6671])
    R={f for f in R if f in H.freeinp}
    # find fragile guards under random combo of R (2 iterations to grow)
    for _ in range(2):
        random.seed(11)
        dv=tangent({f:random.randrange(1,p) for f in R},val)
        vd=[Dm(val[i],dv[i]) for i in range(H.NVARS)]
        ns={'v':vd,'__builtins__':{}}
        frag=[]
        for i in range(len(H.eqcode)):
            if i in Fset: continue
            r=eval(H.eqcode[i],ns)
            if isinstance(r,Dm) and r.b%p!=0: frag.append(i)
        newf=set()
        for i in frag:
            for v in H.eqvars[i]:
                if v in H.freeinp: newf.add(v)
        if newf<=R: break
        R|=newf
    frees=sorted(R)
    nfrag=len(frag)
    # build Jacobian
    M=[[0]*len(frees) for _ in range(nfrag+2)]
    for j,f in enumerate(frees):
        dv=tangent({f:1},val)
        vd=[Dm(val[i],dv[i]) for i in range(H.NVARS)]
        ns={'v':vd,'__builtins__':{}}
        for r,i in enumerate(frag):
            res=eval(H.eqcode[i],ns)
            if isinstance(res,Dm) and res.b: M[r][j]=res.b%p
        M[nfrag][j]=dv[35389]%p; M[nfrag+1][j]=dv[6671]%p
    S0=val[35389]%p; T0=val[6671]%p
    rhs=[0]*nfrag+[(-S0)%p,(-T0)%p]
    # solve M x = rhs (min via echelon, free vars=0)
    nrow=nfrag+2; ncol=len(frees)
    A=[M[r][:]+[rhs[r]] for r in range(nrow)]
    pivinfo=[]; r=0
    for c in range(ncol):
        piv=-1
        for rr in range(r,nrow):
            if A[rr][c]%p!=0: piv=rr;break
        if piv<0: continue
        A[r],A[piv]=A[piv],A[r]
        iv=inv(A[r][c]); A[r]=[(x*iv)%p for x in A[r]]
        for rr in range(nrow):
            if rr!=r and A[rr][c]%p!=0:
                f=A[rr][c]; A[rr]=[(A[rr][k]-f*A[r][k])%p for k in range(ncol+1)]
        pivinfo.append((r,c)); r+=1
        if r==nrow: break
    # inconsistency?
    for rr in range(nrow):
        if all(A[rr][c]%p==0 for c in range(ncol)) and A[rr][ncol]%p!=0:
            return None, frees, nfrag, "INCONSISTENT"
    x=[0]*ncol
    for (rr,c) in pivinfo:
        x[c]=A[rr][ncol]%p
    return x, frees, nfrag, "ok"
# Newton loop
d0=H.loadd('best/new_instance_partial_39013.json')
val=[d0.get(v,0) for v in range(H.NVARS)]
for v in range(H.NVARS): H.val[v]=val[v]
H.forward(); val=list(H.val)
for it in range(8):
    Fset=set(H.fails())
    Score=[i for i in Fset]
    t0=time.time()
    x,frees,nfrag,status=build_and_solve(val,Fset)
    if status=="INCONSISTENT":
        print(f"iter {it}: fails={len(Fset)} HEAL INCONSISTENT (nfrag={nfrag})")
        break
    # apply
    for j,f in enumerate(frees):
        val[f]=(val[f]+x[j])
    for v in range(H.NVARS): H.val[v]=val[v]
    H.forward(); val=list(H.val)
    F=H.fails()
    print(f"iter {it}: applied heal (nfrag={nfrag},nfree={len(frees)}); now fails={len(F)}  ({time.time()-t0:.0f}s)")
    if len(F)==0: 
        print("SOLVED core!"); break
