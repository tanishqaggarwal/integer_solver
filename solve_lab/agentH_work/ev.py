"""Fast forward evaluator over the syntactic gate DAG + exact equation evaluation."""
import model, pickle, os, re, ast, json, sys
from collections import defaultdict
HERE=os.path.dirname(os.path.abspath(__file__))
VAR_RE=re.compile(r'x_(\d+)')
d=model.get(); atom_src=d['atom_src']; atom_vars=d['atom_vars']; eq_terms=d['eq_terms']
F=pickle.load(open(os.path.join(HERE,'fwd2.pkl'),'rb'))
tgt=F['tgt']; definer=F['definer']; order=F['order']; free0=F['free0']; checks=set(F['checks'])
NV=38748

# compile definition expressions: x_t = t - atom  (atom = x_t - rest, so rest = x_t - atom)
def_expr={}
for v in order:
    a=definer[v]
    s=atom_src[a]
    t=ast.parse(s,mode='eval').body
    assert isinstance(t,ast.BinOp) and isinstance(t.op,ast.Sub) and isinstance(t.left,ast.Name)
    rhs=ast.unparse(t.right)
    def_expr[v]=compile(VAR_RE.sub(r'v[\1]', rhs),'<d>','eval')
DEFSEQ=[(v,def_expr[v]) for v in order]
CHECKCODE={a:compile(VAR_RE.sub(r'v[\1]',atom_src[a]),'<c>','eval') for a in checks}

def forward(freevals):
    v=[0]*NV
    for u,val in freevals.items(): v[u]=val
    ns={'v':v,'__builtins__':{}}
    for t,c in DEFSEQ: v[t]=eval(c,ns)
    return v

def eqvals(v):
    ns={'v':v,'__builtins__':{}}
    av={a:eval(c,ns) for a,c in CHECKCODE.items()}
    out=[]
    for m,sq,tl in eq_terms:
        s=0
        for c,a in tl:
            if a in av: s+=c*av[a]
        out.append(m*(s*s if sq else s))
    return out,av

if __name__=='__main__':
    import time
    t0=time.time()
    v=forward({})
    print('forward from all-zero free inputs: %.1fs'%(time.time()-t0))
    ev,av=eqvals(v)
    fails=[i for i,x in enumerate(ev) if x!=0]
    print('failing:',len(fails),'score',len(ev)-len(fails))
    nz=[a for a,x in av.items() if x!=0]
    print('nonzero check atoms:',len(nz))
    json.dump({('x_%d'%i):v[i] for i in range(NV) if v[i]!=0}, open(os.path.join(HERE,'allzero_fwd.json'),'w'))
