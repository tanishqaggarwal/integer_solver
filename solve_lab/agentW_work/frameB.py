"""Parametrised frame: fwd2 orientation with a chosen DETACH set turned back into free inputs."""
import model, pickle, os, re, ast, json
from collections import defaultdict
HERE=os.path.dirname(os.path.abspath(__file__))
VAR_RE=re.compile(r'x_(\d+)')
d=model.get(); atom_src=d['atom_src']; atom_vars=d['atom_vars']; eq_terms=d['eq_terms']
F=pickle.load(open(os.path.join(HERE,'fwd2.pkl'),'rb'))
NV=38748; NEQ=len(eq_terms)
ACODE={a:compile(VAR_RE.sub(r'v[\1]',atom_src[a]),'<a>','eval') for a in range(len(atom_src))}
DEFEXPR={}
for v in F['order']:
    a=F['definer'][v]
    t=ast.parse(atom_src[a],mode='eval').body
    DEFEXPR[v]=compile(VAR_RE.sub(r'v[\1]',ast.unparse(t.right)),'<d>','eval')

class Frame:
    def __init__(s,detach=()):
        s.det=set(detach)
        s.order=[v for v in F['order'] if v not in s.det]
        s.free=sorted(set(F['free0'])|s.det)
        s.definer=F['definer']
        used=set(s.definer[v] for v in s.order)
        s.checks=[a for a in range(len(atom_src)) if a not in used]
        s.fidx={v:i for i,v in enumerate(s.free)}
        sup=[0]*NV
        for v in s.free: sup[v]=1<<s.fidx[v]
        for v in s.order:
            a=s.definer[v]; t=0
            for u in atom_vars[a]:
                if u!=v: t|=sup[u]
            sup[v]=t
        s.sup=sup
        s.csup={a:0 for a in s.checks}
        for a in s.checks:
            t=0
            for u in atom_vars[a]: t|=sup[u]
            s.csup[a]=t
        s.pos={v:i for i,v in enumerate(s.order)}
        s.desc=defaultdict(list)
        for v in s.order:
            t=sup[v]
            while t:
                b=t&-t; s.desc[s.free[b.bit_length()-1]].append(v); t^=b
        for u in s.desc: s.desc[u].sort(key=lambda v:s.pos[v])
        s.chk=defaultdict(list)
        for a in s.checks:
            t=s.csup[a]
            while t:
                b=t&-t; s.chk[s.free[b.bit_length()-1]].append(a); t^=b
        s.eq_of=defaultdict(list)
        for i,(m,sq,tl) in enumerate(eq_terms):
            for c,a in tl: s.eq_of[a].append(i)
        s.SUPV={a:[s.free[b] for b in _bits(s.csup[a])] for a in s.checks}
    def forward(s,fv):
        v=[0]*NV
        for u,val in fv.items(): v[u]=val
        ns={'v':v,'__builtins__':{}}
        for t in s.order: v[t]=eval(DEFEXPR[t],ns)
        return v

def _bits(x):
    o=[]
    while x:
        b=x&-x; o.append(b.bit_length()-1); x^=b
    return o

class State:
    __slots__=('f','fv','v','av','eq','fails','ns')
    def __init__(s,fr,fv):
        s.f=fr; s.fv=dict(fv); s.v=fr.forward(s.fv); s.ns={'v':s.v,'__builtins__':{}}
        s.av={a:eval(ACODE[a],s.ns) for a in fr.checks}
        s.eq=[0]*NEQ
        for i,(m,sq,tl) in enumerate(eq_terms):
            t=0
            for c,a in tl:
                x=s.av.get(a)
                if x: t+=c*x
            s.eq[i]=m*(t*t if sq else t)
        s.fails=set(i for i in range(NEQ) if s.eq[i])
    def clone(s):
        o=State.__new__(State); o.f=s.f; o.fv=dict(s.fv); o.v=s.v[:]
        o.ns={'v':o.v,'__builtins__':{}}; o.av=dict(s.av); o.eq=s.eq[:]; o.fails=set(s.fails)
        return o
    def set_free(s,ch):
        fr=s.f; aff=set(); ck=set()
        for u,val in ch.items():
            if s.v[u]==val and s.fv.get(u,0)==val: continue
            s.fv[u]=val; s.v[u]=val; aff.update(fr.desc[u]); ck.update(fr.chk[u])
        for v in sorted(aff,key=lambda v:fr.pos[v]): s.v[v]=eval(DEFEXPR[v],s.ns)
        eqs=set()
        for a in ck:
            nv=eval(ACODE[a],s.ns)
            if nv!=s.av[a]: s.av[a]=nv; eqs.update(fr.eq_of[a])
        for i in eqs:
            m,sq,tl=eq_terms[i]; t=0
            for c,a in tl:
                x=s.av.get(a)
                if x: t+=c*x
            val=m*(t*t if sq else t); s.eq[i]=val
            if val: s.fails.add(i)
            else: s.fails.discard(i)
        return s
    def score(s): return NEQ-len(s.fails)
    def nz(s): return [a for a,x in s.av.items() if x]

if __name__=='__main__':
    import time
    DET=[642,28730,29854,31864]
    t0=time.time(); fr=Frame(DET); print('frame built %.1fs free=%d checks=%d'%(time.time()-t0,len(fr.free),len(fr.checks)))
    W=json.load(open('../best/new_instance_partial_39026.json'))
    v=[0]*NV
    for k,val in W.items(): v[int(k[2:]) if k.startswith('x_') else int(k)]=int(val)
    fv={u:v[u] for u in fr.free if v[u]!=0}
    st=State(fr,fv)
    print('witness in detached frame: score',st.score(),'nz',sorted(st.nz()),'failing',sorted(st.fails))
    dif=[i for i in range(NV) if st.v[i]!=v[i]]
    print('vars differing from witness:',len(dif),dif[:10])
    pickle.dump(DET,open('frameB_det.pkl','wb'))
