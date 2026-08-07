#!/usr/bin/env python3
"""AUDIT T26 -- agent N's detach exhaustion, the largest 'reported' row in my own ledger.

N: 'make(D) gives detached pool members their WITNESS values.  Only 4 of the 65 pool variables
have a witness value different from their gate value -- {642,28730,29854,31864}.  Detaching any
of the other 61 is a literal no-op.  So the whole 2^65 lattice has 16 distinct states.'

That reduction IS the claim.  Test it in F's certified-faithful parse (audit T2), not N's.
KEY IDENTITY: a pool variable v is defined by an atom (v - RHS).  witness(v) != gate(v) at the
deliverable  <=>  that atom is NONZERO at the deliverable.  I already know F's parse has exactly
7 nonzero atoms there, so this is directly checkable."""
import os,sys,json,pickle,collections,re
T=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.abspath(os.path.join(T,'..'))
F=os.path.join(LAB,'agentF_work'); sys.path.insert(0,F); sys.path.insert(0,LAB)
from fwd import compile_node
from circ2 import vars_of
d=pickle.load(open(os.path.join(F,'circ4.pkl'),'rb'))
atoms=d['atoms']; names=list(atoms)
POOL=json.load(open(os.path.join(LAB,'agentN_work','pool.json')))
WIT=[642,28730,29854,31864]
print('pool size: %d   N\'s witness set: %s'%(len(POOL),WIT))
NV=38748
B=json.load(open(os.path.join(LAB,'best','new_instance_partial_39026.json')))
v=[0]*NV
for k,val in B.items(): v[int(k[2:])]=int(val)
prog=compile('r[:]=['+','.join(compile_node(atoms[a]) for a in names)+']','<at>','exec')
r=[0]*len(names); exec(prog,{'v':v,'r':r,'__builtins__':{}})
nz={names[i] for i in range(len(names)) if r[i]}
print('nonzero atoms at the deliverable (F parse): %d'%len(nz))
# for each pool var, find its DEFINING atom: shape (xV - RHS) or (RHS - xV) with V the defined var
defatom={}
for i,a in enumerate(names):
    s=a.replace(' ','')
    m=re.match(r'^\(x(\d+)-(.+)\)$',s)
    if m:
        vv=int(m.group(1))
        if vv in set(POOL) and vv not in defatom: defatom[vv]=(a,r[i])
missing=[u for u in POOL if u not in defatom]
print('pool vars with an identifiable defining atom (xV - RHS): %d ; without: %d'%(len(defatom),len(missing)))
diff=[u for u in defatom if defatom[u][1]!=0]
print('\n== pool vars whose defining atom is NONZERO at the deliverable (witness != gate) ==')
for u in sorted(diff): print('   x%-6d  %s'%(u,defatom[u][0][:70]))
print('   count: %d   equals N\'s witness set exactly? %s'%(len(diff),sorted(diff)==sorted(WIT)))
print('   the other %d pool vars have their defining atom ZERO -> re-attaching is a no-op AT THE WITNESS STATE'%(len(defatom)-len(diff)))
# --- the subtlety: does any non-WIT pool var DEPEND on a WIT member?  If so its gate value
#     changes when a WIT member is re-attached, and "no-op" needs to hold at 16 states, not one.
print('\n== do the other 61 depend on the witness set? ==')
dep=collections.defaultdict(set)
for u,(a,val) in defatom.items():
    dep[u]=vars_of(atoms[a])-{u}
direct=[u for u in defatom if u not in WIT and dep[u]&set(WIT)]
print('   directly referencing a WIT variable: %d  %s'%(len(direct),direct[:8]))
# transitive within the pool
poolset=set(POOL); reach={u:set(dep[u]&poolset) for u in defatom}
for _ in range(8):
    for u in reach:
        add=set()
        for w in list(reach[u]):
            add|=reach.get(w,set())
        reach[u]|=add
trans=[u for u in defatom if u not in WIT and reach[u]&set(WIT)]
print('   transitively (within the pool) reaching a WIT variable: %d  %s'%(len(trans),trans[:8]))
if not trans:
    print('   -> NONE.  The 61 gate values cannot change with D, so the no-op holds at all 16')
    print('      states, not just at the witness state.  The reduction to 16 is a PROOF.')
else:
    print('   -> some do; the no-op must then be checked at all 16 states, not just the witness.')

print('\n== CLOSING THE GAP: the check above was POOL-INTERNAL.  A non-WIT pool var could reach')
print('   a WIT var through variables OUTSIDE the pool.  Redo over the whole definition DAG. ==')
alldef={}
for i,a in enumerate(names):
    s=a.replace(' ','')
    m=re.match(r'^\(x(\d+)-(.+)\)$',s)
    if m:
        vv=int(m.group(1))
        if vv not in alldef: alldef[vv]=vars_of(atoms[a])-{vv}
print('   variables with a definition atom of shape (xV - RHS): %d'%len(alldef))
WITS=set(WIT)
def reaches_wit(start):
    seen=set(); st=list(start)
    while st:
        u=st.pop()
        if u in seen: continue
        seen.add(u)
        if u in WITS: return True
        # a DETACHED pool var is free -> its definition is not used; but for the worst case we
        # follow every definition edge regardless, which can only OVER-report dependence
        for w in alldef.get(u,()): 
            if w not in seen: st.append(w)
    return False
bad=[]
for u in defatom:
    if u in WITS: continue
    if reaches_wit(dep[u]): bad.append(u)
print('   of the 61, reaching a WIT variable anywhere in the full DAG: %d  %s'%(len(bad),bad[:10]))
if not bad:
    print('   -> STILL NONE, over the whole instance.  The 61 gate values are independent of D.')
    print('      make(D) depends only on D & {642,28730,29854,31864}: the 2^65 lattice has')
    print('      EXACTLY 16 states, by proof.  The 16 measured signatures are complete.')
else:
    print('   -> the reduction needs the no-op verified at all 16 states, not just the witness.')
