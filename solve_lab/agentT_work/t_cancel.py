#!/usr/bin/env python3
"""AUDIT T15 -- agent L's cancellation result, the basis on which M and others redirected.

L: with BYTE-IDENTICAL support (the deliverable's own 4 atoms) the deliverable costs 7 failing
equations and L's build2 costs 13; a 2-atom variant costs 11.  L attributes the whole 6-equation
gap to 12 cofactor/handle variables the deliverable sets nonzero and L's constructor leaves at 0:
  x105 x1329 x3387 x5081 x5676 x9413 x10903 x11436 x14393 x14768 x17325 x22820

Tested from the DELIVERABLE side so the result does not depend on L's constructor (which L itself
flags as carrying an un-converged divisibility repair that adds unrelated nonzero atoms):
  zero those 12 in the deliverable and measure BOTH the score AND the nonzero-atom support.
If support is unchanged and the score moves 7 -> 13, 'cost is a value property' is established.
Atom support is read with agent F's decomposition, which I certified faithful in audit T2."""
import os,sys,json,pickle,collections,subprocess
T=os.path.dirname(os.path.abspath(__file__))
F=os.path.abspath(os.path.join(T,'..','agentF_work')); sys.path.insert(0,F)
sys.path.insert(0,os.path.abspath(os.path.join(T,'..')))
from fwd import compile_node
import checker as CK
d=pickle.load(open(os.path.join(F,'circ4.pkl'),'rb'))
atoms=d['atoms']; names=list(atoms)
prog=compile('r[:]=['+','.join(compile_node(atoms[a]) for a in names)+']','<at>','exec')
NV=38748
BASE=json.load(open(os.path.join(T,'..','best','new_instance_partial_39026.json')))
V=[0]*NV
for k,val in BASE.items(): V[int(k[2:])]=int(val)
TWELVE=[105,1329,3387,5081,5676,9413,10903,11436,14393,14768,17325,22820]
eqs=CK.load_equations()
print('equations loaded: %d'%len(eqs),flush=True)
# which of the 12 are FREE variables (never on the left of a definition) -- use F's definer if any
import circ2
defined=set()
for a in names:
    pass
def score_and_support(v):
    r=[0]*len(names); exec(prog,{'v':v,'r':r,'__builtins__':{}})
    sup=sorted(names[i] for i in range(len(names)) if r[i])
    return sup
def checker_score(v,tag):
    asg={'x_%d'%i:str(v[i]) for i in range(NV) if v[i]}
    fn=os.path.join(T,'_tmp_%s.json'%tag); json.dump(asg,open(fn,'w'))
    out=subprocess.run([sys.executable,os.path.join(T,'..','checker.py'),fn],capture_output=True,text=True).stdout
    os.remove(fn)
    return out.strip().splitlines()[-1] if out.strip() else '(no output)'
print('\n-- deliverable as given --',flush=True)
sup0=score_and_support(V)
print('   nonzero atoms: %d  %s'%(len(sup0),sup0),flush=True)
print('   checker: %s'%checker_score(V,'d0'),flush=True)
print('\n-- deliverable with L\'s 12 cofactor vars ZEROED --',flush=True)
V2=list(V)
for u in TWELVE: V2[u]=0
sup2=score_and_support(V2)
print('   nonzero atoms: %d  %s'%(len(sup2),sup2),flush=True)
print('   support identical to deliverable? %s'%(sup2==sup0),flush=True)
if sup2!=sup0:
    print('   added: %s'%sorted(set(sup2)-set(sup0)),flush=True)
    print('   removed: %s'%sorted(set(sup0)-set(sup2)),flush=True)
print('   checker: %s'%checker_score(V2,'d12'),flush=True)
print('\n-- one at a time: zero each of the 12 alone --',flush=True)
for u in TWELVE:
    V3=list(V); V3[u]=0
    s=score_and_support(V3)
    print('   x%-6d value %-14s -> nonzero atoms %2d (same support: %-5s)  %s'%(
        u,str(V[u])[:14],len(s),s==sup0,checker_score(V3,'d%d'%u)),flush=True)
