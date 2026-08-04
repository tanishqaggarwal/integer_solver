import heal_harness as H, pickle, time, random
p=H.p
AD=pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/adcache.pkl','rb'))
gate_parsed=AD['gate_parsed']; val=AD['val']; gS=AD['gS']; gT=AD['gT']
for v in range(H.NVARS): H.val[v]=val[v]
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
def tangent(seedmap):
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
# iter-1 closure: movers -> fragile guards -> frees
R0=sorted(f for f in H.freeinp if gS[f]%p or gT[f]%p)
random.seed(7)
dv=tangent({f:random.randrange(1,p) for f in R0})
vdual=[Dm(val[i],dv[i]) for i in range(H.NVARS)]
ns={'v':vdual,'__builtins__':{}}
Fset=set(H.fails())
frag=[i for i in range(len(H.eqcode)) if i not in Fset and isinstance(eval(H.eqcode[i],ns),Dm) and eval(H.eqcode[i],ns).b%p!=0]
frees=set(R0)
for i in frag:
    for v in H.eqvars[i]:
        if v in H.freeinp: frees.add(v)
frees=sorted(frees)
print(f"closure: {len(frees)} frees, {len(frag)} fragile guards")
fidx={f:j for j,f in enumerate(frees)}
# Build Jacobian rows: for each free (column) compute tangent and each guard/S/T derivative
# We'll build dense matrix M with rows=guards+[S,T], cols=frees
nrow=len(frag)+2; ncol=len(frees)
M=[[0]*ncol for _ in range(nrow)]
t0=time.time()
for j,f in enumerate(frees):
    dv=tangent({f:1})
    vd=[Dm(val[i],dv[i]) for i in range(H.NVARS)]
    ns={'v':vd,'__builtins__':{}}
    for r,i in enumerate(frag):
        res=eval(H.eqcode[i],ns)
        if isinstance(res,Dm) and res.b: M[r][j]=res.b%p
    # S,T
    M[len(frag)][j]=dv[35389]%p
    M[len(frag)+1][j]=dv[6671]%p
print(f"built Jacobian {nrow}x{ncol} in {time.time()-t0:.1f}s")
# RHS: guards target 0, S target -S, T target -T
S0=val[35389]%p; T0=val[6671]%p
rhs=[0]*len(frag)+[(-S0)%p,(-T0)%p]
pickle.dump({'M':M,'rhs':rhs,'frees':frees,'nfrag':len(frag)},open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/jacmat.pkl','wb'))
print("saved. S0=",S0!=0,"T0=",T0!=0)
