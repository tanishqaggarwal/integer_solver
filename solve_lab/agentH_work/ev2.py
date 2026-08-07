"""Evaluator on the extended (fwd4) orientation.  Generic unit-coefficient definitions."""
import model, pickle, os, re, json
from collections import defaultdict
HERE=os.path.dirname(os.path.abspath(__file__))
VAR_RE=re.compile(r'x_(\d+)')
d=model.get(); atom_src=d['atom_src']; atom_vars=d['atom_vars']; eq_terms=d['eq_terms']
polys=pickle.load(open(os.path.join(HERE,'polys.pkl'),'rb'))
F4=pickle.load(open(os.path.join(HERE,'fwd4.pkl'),'rb'))
definer=F4['definer']; order=F4['order']; free=F4['free']; checks=F4['checks']
NV=38748
ACODE={a:compile(VAR_RE.sub(r'v[\1]',atom_src[a]),'<a>','eval') for a in range(len(atom_src))}
COEF={}
for v in order:
    a=definer[v]; COEF[v]=polys[a][(v,)]
DEFSEQ=[(v,definer[v],COEF[v]) for v in order]
CHECKS=list(checks)
def forward(fv):
    v=[0]*NV
    for u,val in fv.items(): v[u]=val
    ns={'v':v,'__builtins__':{}}
    for t,a,c in DEFSEQ:
        v[t]=0
        v[t]=-eval(ACODE[a],ns)*c
    return v
def eqvals(v):
    ns={'v':v,'__builtins__':{}}
    av={a:eval(ACODE[a],ns) for a in CHECKS}
    out=[]
    for m,sq,tl in eq_terms:
        s=0
        for c,a in tl:
            x=av.get(a)
            if x: s+=c*x
        out.append(m*(s*s if sq else s))
    return out,av
if __name__=='__main__':
    import time
    t0=time.time(); v=forward({}); e,av=eqvals(v)
    f=[i for i,x in enumerate(e) if x!=0]
    print('fwd4 all-zero: score',39033-len(f),'nz',sorted(a for a,x in av.items() if x),'%.1fs'%(time.time()-t0))
    print('free inputs',len(free),'checks',len(CHECKS))
