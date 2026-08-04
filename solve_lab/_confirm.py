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
# 3-iteration closure with seed 999
R={f for f in (H.anc[35389]|H.anc[6671]) if f in H.freeinp}
for it in range(3):
    random.seed(999+it)
    dv=tangent({f:random.randrange(1,p) for f in R},val)
    vd=[Dm(val[i],dv[i]) for i in range(H.NVARS)]; ns={'v':vd,'__builtins__':{}}
    frag=[i for i in range(len(H.eqcode)) if i not in Fset and isinstance(eval(H.eqcode[i],ns),Dm) and eval(H.eqcode[i],ns).b%p!=0]
    newf=set()
    for i in frag:
        for v in H.eqvars[i]:
            if v in H.freeinp: newf.add(v)
    if newf<=R: break
    R|=newf
frees=sorted(R); nfrag=len(frag)
print(f"seed999 closure {len(frees)} frees, {nfrag} guards")
# test ONLY: is dS in rowspace(guards)? is dT? Solve by augmenting guards with S row and T row separately.
def build_cols():
    cols=[]
    for f in frees:
        dv=tangent({f:1},val)
        vd=[Dm(val[i],dv[i]) for i in range(H.NVARS)]; ns={'v':vd,'__builtins__':{}}
        col=[0]*(nfrag+2)
        for r,i in enumerate(frag):
            res=eval(H.eqcode[i],ns)
            if isinstance(res,Dm) and res.b: col[r]=res.b%p
        col[nfrag]=dv[35389]%p; col[nfrag+1]=dv[6671]%p
        cols.append(col)
    return cols
cols=build_cols()
ncol=len(frees)
def rankof(rowids):
    B=[[cols[j][r] for j in range(ncol)] for r in rowids]
    rr=0
    for c in range(ncol):
        piv=-1
        for k in range(rr,len(B)):
            if B[k][c]%p!=0: piv=k;break
        if piv<0: continue
        B[rr],B[piv]=B[piv],B[rr]
        ivv=inv(B[rr][c]); B[rr]=[(x*ivv)%p for x in B[rr]]
        for k in range(len(B)):
            if k!=rr and B[k][c]%p!=0:
                f=B[k][c]; B[k]=[(B[k][j]-f*B[rr][j])%p for j in range(ncol)]
        rr+=1
        if rr==len(B): break
    return rr
g=list(range(nfrag))
rg=rankof(g); rgS=rankof(g+[nfrag]); rgT=rankof(g+[nfrag+1]); rgST=rankof(g+[nfrag,nfrag+1])
print(f"rank(guards)={rg}")
print(f"rank(guards+S)={rgS} -> S {'LOCKED' if rgS==rg else 'movable'}")
print(f"rank(guards+T)={rgT} -> T {'LOCKED' if rgT==rg else 'movable'}")
print(f"rank(guards+S+T)={rgST} -> core DOF relative to guards = {rgST-rg}")
