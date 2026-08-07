"""THE decomposition: connected components of the free-input hypergraph induced by equations."""
import model, pickle, os, json, time
from collections import Counter, defaultdict
HERE=os.path.dirname(os.path.abspath(__file__))
d=model.get(); eq_terms=d['eq_terms']
S=pickle.load(open(os.path.join(HERE,'support.pkl'),'rb'))
esup=S['esup']; fidx=S['fidx']; csup=S['csup']
NF=len(fidx)
inv={i:v for v,i in fidx.items()}
def bits(x):
    out=[]
    while x:
        b=x&-x; out.append(b.bit_length()-1); x^=b
    return out
class DSU:
    def __init__(s,n): s.p=list(range(n))
    def f(s,x):
        while s.p[x]!=x: s.p[x]=s.p[s.p[x]]; x=s.p[x]
        return x
    def u(s,a,b):
        a,b=s.f(a),s.f(b)
        if a!=b: s.p[a]=b
dsu=DSU(NF)
eqbits=[]
for i,e in enumerate(esup):
    bs=bits(e); eqbits.append(bs)
    for x in bs[1:]: dsu.u(bs[0],x)
comp=defaultdict(list)
touched=set()
for bs in eqbits: touched.update(bs)
for x in touched: comp[dsu.f(x)].append(x)
sz=sorted((len(v) for v in comp.values()),reverse=True)
print('free inputs total', NF, 'appearing in some non-trivial eq:', len(touched))
print('components:', len(sz), 'sizes top 30:', sz[:30])
print('size hist:', Counter(sz).most_common(15))
# per component: how many equations
compof={}
for r,vs in comp.items():
    for v in vs: compof[v]=r
eqcomp=defaultdict(list)
trivial=0
for i,bs in enumerate(eqbits):
    if not bs: trivial+=1; continue
    eqcomp[compof[bs[0]]].append(i)
print('equations with empty support (auto-satisfied by forward eval):', trivial)
info=sorted(((len(comp[r]), len(eqcomp[r]), r) for r in comp), reverse=True)
print('top components (nfree, neq):', [(a,b) for a,b,_ in info[:20]])
print('total eqs in components:', sum(b for a,b,_ in info))
# where do the 7 failing equations live?
FAIL=[12231,12270,12350,14584,18673,22044,29125]
for f in FAIL:
    bs=eqbits[f]
    r=compof[bs[0]] if bs else None
    print('fail eq',f,'support',len(bs),'comp size',len(comp[r]) if r is not None else 0,
          'comp eqs',len(eqcomp[r]) if r is not None else 0)
json.dump({'sizes':sz,'ncomp':len(sz)},open(os.path.join(HERE,'decomp_sizes.json'),'w'))
pickle.dump({'comp':dict(comp),'eqcomp':dict(eqcomp),'compof':compof,'eqbits':eqbits,'inv':inv},
            open(os.path.join(HERE,'decomp.pkl'),'wb'))
