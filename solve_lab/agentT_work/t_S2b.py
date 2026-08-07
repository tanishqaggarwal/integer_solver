#!/usr/bin/env python3
"""AUDIT T24b -- is L's |S|=2 closure true OF THE INSTANCE, not just of L's 9,032-atom engine?
Evaluate all 39,033 atoms of F's certified-faithful decomposition at the reproduced assignment
and check that the nonzero set is exactly the two target congruences, and that their equation
footprint is exactly the 15 failing lines checker.py reports."""
import os,sys,json,pickle,collections
T=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.abspath(os.path.join(T,'..'))
F=os.path.join(LAB,'agentF_work'); sys.path.insert(0,F); sys.path.insert(0,LAB)
from fwd import compile_node
import checker as CK
d=pickle.load(open(os.path.join(F,'circ4.pkl'),'rb'))
atoms=d['atoms']; eqrows=d['eqrows']; names=list(atoms); idx={a:i for i,a in enumerate(names)}
a2e=collections.defaultdict(set)
for e,row in enumerate(eqrows):
    for k,a in row: a2e[idx[a]].add(e)
prog=compile('r[:]=['+','.join(compile_node(atoms[a]) for a in names)+']','<at>','exec')
NV=38748
codes,varsets=CK.load_equations()
def load(fn):
    v=[0]*NV
    for k,val in json.load(open(fn)).items(): v[int(k[2:])]=int(val)
    return v
for tag,fn in [('L |S|=2 (reproduced by T)',os.path.join(T,'t_S2_assign.json')),
               ('L |S|=1 assign_L1.json',os.path.join(LAB,'agentL_work','assign_L1.json'))]:
    v=load(fn)
    r=[0]*len(names); exec(prog,{'v':v,'r':r,'__builtins__':{}})
    nz=[i for i in range(len(names)) if r[i]]
    fails=CK.evaluate_all(codes,v)
    foot=set()
    for i in nz: foot|=a2e[i]
    print('\n== %s =='%tag)
    print('   nonzero atoms in F\'s parse : %d'%len(nz))
    for i in nz: print('        %s'%names[i][:80])
    print('   checker failing            : %d  %s'%(len(fails),fails))
    print('   equation footprint of those atoms : %d'%len(foot))
    print('   footprint == failing set ? %s'%(foot==set(fails)))
    if foot!=set(fails):
        print('     in footprint not failing (cancellation): %s'%sorted(foot-set(fails)))
        print('     failing not in footprint (UNEXPLAINED) : %s'%sorted(set(fails)-foot))
