import heal_harness as H, pickle, time
p=H.p
AD=pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/adcache.pkl','rb'))
SC=pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/satcache.pkl','rb'))
gate_parsed=AD['gate_parsed']; val=AD['val']; gS=AD['gS']; gT=AD['gT']
sat=SC['sat']; var2eq=SC['var2eq']
satset=set(sat)
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
def combined_tangent(seedmap):
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
def fragile_under(R):
    import random; random.seed(7)
    dv=combined_tangent({f:random.randrange(1,p) for f in R})
    vdual=[Dm(val[i],dv[i]) for i in range(H.NVARS)]
    ns={'v':vdual,'__builtins__':{}}
    frag=[]
    for i in sat:
        r=eval(H.eqcode[i],ns)
        if isinstance(r,Dm) and r.b%p!=0: frag.append(i)
    return frag
R=set(f for f in H.freeinp if gS[f]%p or gT[f]%p)
for it in range(6):
    t0=time.time()
    frag=fragile_under(R)
    # expand R with free vars in fragile guards
    newfrees=set()
    for i in frag:
        for v in H.eqvars[i]:
            if v in H.freeinp: newfrees.add(v)
    grow=newfrees-R
    print(f"iter {it}: |R|={len(R)} fragile={len(frag)} new_frees={len(grow)} ({time.time()-t0:.1f}s)")
    if not grow: break
    R|=grow
print(f"\nCLOSURE: frees={len(R)} fragile_guards={len(frag)}")
pickle.dump({'R':sorted(R),'frag':sorted(frag)}, open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/closurecache.pkl','wb'))
