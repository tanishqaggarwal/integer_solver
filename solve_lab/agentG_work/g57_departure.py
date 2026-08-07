"""Coset decoding, step 3: evaluate a departure supported on a chosen small set of
unknowns.  Everything else is held at the all-linear-forced solution x0.  Each equation
becomes an exact polynomial in the departure variables; count how many can be zeroed."""
import os, sys, pickle, itertools, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G
from gsym2 import L, ad, P
D=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/coset_model.pkl','rb'))
Lin=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/coset_lin.pkl','rb'))
NB=D['NB']; lin=D['lin']; non=D['non']; n=len(NB); x0=Lin['x0']
ix={u:i for i,u in enumerate(NB)}
SUP=[ix[int(u)] for u in sys.argv[1].split(',')]
print('departure support:',[NB[c] for c in SUP])
k=len(SUP); sidx={c:j for j,c in enumerate(SUP)}
def sub(f):
    """polynomial in the departure variables only (others fixed at x0)"""
    out={}
    for m,c in f.items():
        t=c%P; e=[0]*k
        for col,ee in m:
            if col in sidx: e[sidx[col]]+=ee
            else: t=t*pow(x0[col],ee,P)%P
        if t: 
            key=tuple(e); out[key]=(out.get(key,0)+t)%P
    return {m:c for m,c in out.items() if c}
polys=[]
for i,f in lin+ [ (i,(f if not isinstance(f,int) else {():f})) for i,f in non]:
    g=sub(f if not isinstance(f,int) else {():f})
    polys.append((i,g))
const_bad=[i for i,g in polys if len(g)==1 and () in g]
varying=[(i,g) for i,g in polys if any(any(e) for e in g)]
zero=[i for i,g in polys if not g]
print('equations: %d total  |  identically zero %d  |  nonzero CONSTANT %d  |  varying %d'
      %(len(polys),len(zero),len(const_bad),len(varying)))
print('nonzero constants (unfixable by this departure):',const_bad[:30])
print('varying equations:',[i for i,_ in varying][:40])
for i,g in varying[:40]:
    terms=' + '.join('%s*%s'%(('%d'%c if c<10**10 else 'C'),('*'.join('d%d^%d'%(NB[SUP[j]],e) for j,e in enumerate(m) if e) or '1')) for m,c in sorted(g.items(),key=lambda kv:-sum(kv[0])))
    print('   eq%-6d %s'%(i,terms[:150]))
pickle.dump({'polys':polys,'SUP':[NB[c] for c in SUP]},open('/home/user/integer_solver/solve_lab/agentG_work/dep_%s.pkl'%sys.argv[1].replace(',','_'),'wb'))
