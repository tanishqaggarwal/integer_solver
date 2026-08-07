#!/usr/bin/env python3
"""AUDIT T15 -- agent L's cancellation result, the basis on which M redirected its search.

L: with BYTE-IDENTICAL support (the deliverable's own 4 atoms) the deliverable costs 7 failing
equations and L's build2 costs 13; a 2-atom variant costs 11.  L attributes the whole 6-equation
gap to 12 cofactor/handle vars the deliverable sets nonzero and L's constructor leaves at 0:
  x105 x1329 x3387 x5081 x5676 x9413 x10903 x11436 x14393 x14768 x17325 x22820

Tested from the DELIVERABLE side, so the result does not depend on L's constructor -- which L
itself flags as carrying an un-converged divisibility repair that adds unrelated nonzero atoms.
Zero the 12 in the deliverable and measure BOTH the exact score AND the nonzero-atom support.
Atom support is read in agent F's decomposition, certified faithful in audit T2."""
import os,sys,json,pickle,time
T=os.path.dirname(os.path.abspath(__file__))
LAB=os.path.abspath(os.path.join(T,'..'))
F=os.path.join(LAB,'agentF_work'); sys.path.insert(0,F); sys.path.insert(0,LAB)
sys.path.insert(0,os.path.join(LAB,'agentE_work'))
from fwd import compile_node
import checker as CK
d=pickle.load(open(os.path.join(F,'circ4.pkl'),'rb'))
atoms=d['atoms']; names=list(atoms)
prog=compile('r[:]=['+','.join(compile_node(atoms[a]) for a in names)+']','<at>','exec')
NV=38748
t0=time.time(); codes,varsets=CK.load_equations(); print('equations: %d (%.0fs)'%(len(codes),time.time()-t0),flush=True)
BASE=json.load(open(os.path.join(LAB,'best','new_instance_partial_39026.json')))
V=[0]*NV
for k,val in BASE.items(): V[int(k[2:])]=int(val)
TWELVE=[105,1329,3387,5081,5676,9413,10903,11436,14393,14768,17325,22820]
import engine as E
print('\nare the 12 FREE variables (no definition) in agent E\'s parse?')
for u in TWELVE:
    print('   x%-6d free=%-5s   deliverable value = %s'%(u,E.definer[u] is None,str(V[u])[:56]))
print('   -> free: %d of 12'%sum(1 for u in TWELVE if E.definer[u] is None),flush=True)
r=[0]*len(names)
def sup_of(v):
    exec(prog,{'v':v,'r':r,'__builtins__':{}})
    return tuple(names[i] for i in range(len(names)) if r[i])
def fails(v): return CK.evaluate_all(codes,v)
sup0=sup_of(V); f0=fails(V)
print('\n-- deliverable as given --',flush=True)
print('   nonzero atoms %d: %s'%(len(sup0),list(sup0)),flush=True)
print('   FAILING %d: %s'%(len(f0),f0),flush=True)
V2=list(V)
for u in TWELVE: V2[u]=0
sup2=sup_of(V2); f2=fails(V2)
print('\n-- deliverable with L\'s 12 cofactor vars ZEROED --',flush=True)
print('   nonzero atoms %d: %s'%(len(sup2),list(sup2)),flush=True)
print('   support IDENTICAL to deliverable? %s'%(sup2==sup0),flush=True)
if sup2!=sup0:
    print('     added  : %s'%sorted(set(sup2)-set(sup0)),flush=True)
    print('     removed: %s'%sorted(set(sup0)-set(sup2)),flush=True)
print('   FAILING %d: %s'%(len(f2),f2),flush=True)
print('   new failures vs deliverable: %s'%sorted(set(f2)-set(f0)),flush=True)
print('\n-- zero each of the 12 alone (which ones actually carry the cancellation?) --',flush=True)
for u in TWELVE:
    V3=list(V); V3[u]=0
    s=sup_of(V3); f=fails(V3)
    print('   x%-6d -> atoms %2d sameSupport=%-5s FAILING %2d   delta %+d'%(
        u,len(s),s==sup0,len(f),len(f)-len(f0)),flush=True)
