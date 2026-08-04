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
d0=H.loadd('best/new_instance_partial_39013.json')
val=[d0.get(v,0) for v in range(H.NVARS)]
for v in range(H.NVARS): H.val[v]=val[v]
H.forward(); val=list(H.val)
Fset=set(H.fails())
# build 2-iter closure
R={f for f in (H.anc[35389]|H.anc[6671]) if f in H.freeinp}
for _ in range(2):
    random.seed(11)
    dv=tangent({f:random.randrange(1,p) for f in R},val)
    vd=[Dm(val[i],dv[i]) for i in range(H.NVARS)]; ns={'v':vd,'__builtins__':{}}
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
frees=sorted(R); nfrag=len(frag)
print(f"closure {len(frees)} frees, {nfrag} guards")
# Build augmented matrix: [jac(frees) | rhs | tagS | tagT]
ncol=len(frees)
rows=[]
t0=time.time()
# to build: iterate frees as columns
Jcols=[]
for j,f in enumerate(frees):
    dv=tangent({f:1},val)
    vd=[Dm(val[i],dv[i]) for i in range(H.NVARS)]; ns={'v':vd,'__builtins__':{}}
    col=[0]*(nfrag+2)
    for r,i in enumerate(frag):
        res=eval(H.eqcode[i],ns)
        if isinstance(res,Dm) and res.b: col[r]=res.b%p
    col[nfrag]=dv[35389]%p; col[nfrag+1]=dv[6671]%p
    Jcols.append(col)
print(f"built in {time.time()-t0:.0f}s")
S0=val[35389]%p; T0=val[6671]%p
# augmented rows: row r -> [Jcols[j][r] for j] + [rhs_r, tagS_r, tagT_r]
A=[]
for r in range(nfrag+2):
    row=[Jcols[j][r] for j in range(ncol)]
    if r<nfrag: row+=[0,0,0]
    elif r==nfrag: row+=[(-S0)%p,1,0]  # S row: target -S, tagS=1
    else: row+=[(-T0)%p,0,1]           # T row
    A.append(row)
W=ncol+3
# Gaussian eliminate on free columns only (0..ncol-1)
nrow=len(A); r=0
for c in range(ncol):
    piv=-1
    for rr in range(r,nrow):
        if A[rr][c]%p!=0: piv=rr;break
    if piv<0: continue
    A[r],A[piv]=A[piv],A[r]
    iv=inv(A[r][c]); A[r]=[(x*iv)%p for x in A[r]]
    for rr in range(nrow):
        if rr!=r and A[rr][c]%p!=0:
            f=A[rr][c]; A[rr]=[(A[rr][k]-f*A[r][k])%p for k in range(W)]
    r+=1
    if r==nrow: break
rank=r
# find rows with all free-cols 0
certs=[]
for rr in range(nrow):
    if all(A[rr][c]%p==0 for c in range(ncol)):
        rhs=A[rr][ncol]%p; lamS=A[rr][ncol+1]%p; lamT=A[rr][ncol+2]%p
        if lamS or lamT:  # a relation involving S/T
            certs.append((rhs,lamS,lamT))
print(f"rank(free part)={rank}")
print(f"conservation certificates (rhs, lamS, lamT):")
for rhs,lamS,lamT in certs:
    # meaning: lamS*dS + lamT*dT is in guard rowspace; and lamS*(-S)+lamT*(-T)=rhs is forced
    # => lamS*S + lamT*T is CONSERVED (= -rhs mod p, i.e. can't be moved from current)
    print(f"   lamS={lamS}, lamT={lamT}, forced target rhs={rhs}  (0 means healable)")
    conserved_val=(lamS*S0+lamT*T0)%p
    print(f"     => {lamS}*S + {lamT}*T is CONSERVED at {conserved_val} (need 0). Movable: {rhs==0}")
if not certs:
    print("   NONE -> S,T both movable (heal consistent)")
