import pickle, re, json, sys, collections
m=pickle.load(open('model3.pkl','rb'))
atoms=m['atoms']; eqt=m['eq_terms']
VAR=re.compile(r'x_(\d+)')
NV=38748
codes=[compile(VAR.sub(r'v[\1]',a),'<a>','eval') for a in atoms]
def load(p):
    d=json.load(open(p)); v=[0]*NV
    for k,val in d.items(): v[int(k[2:]) if k.startswith('x_') else int(k)]=int(val)
    return v
def atomvals(v):
    ns={'v':v,'__builtins__':{}}
    return [eval(c,ns) for c in codes]
def eqvals(av):
    out=[]
    for issq,outer,terms in eqt:
        s=0
        for c,a in terms:
            s+= c*(1 if a<0 else av[a]) if a>=0 else c
        # careful
        out.append(s)
    return out
def eqvals2(av):
    out=[]
    for issq,outer,terms in eqt:
        s=0
        for c,a in terms:
            if a<0: s+=c
            else: s+=c*av[a]
        out.append((s*s*outer if issq else s*outer))
    return out
if __name__=='__main__':
    v=load(sys.argv[1])
    av=atomvals(v)
    nz=[i for i,x in enumerate(av) if x!=0]
    print("nonzero atoms:",len(nz))
    for i in nz: print("  atom",i,atoms[i][:200],"=",av[i])
    ev=eqvals2(av)
    fails=[i for i,x in enumerate(ev) if x!=0]
    print("failing eqs from atom model:",len(fails),fails[:20])
