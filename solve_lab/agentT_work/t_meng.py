#!/usr/bin/env python3
"""AUDIT T27 -- is agent M's incremental engine EXACT, checked outside M's own parse?

M's 2^12 enumeration (4,096 subsets, complete=True, nothing above 39,026) is only as good as its
scorer.  M's six calibration gates include 'incremental == full engine3, 0 vars differing' -- but
that is M checking M.  The cheap independent check is the deliverable's own subset:
drive ieng.tune on H12's witness subset {642,28730,29854,31864}, materialise the assignment M's
engine actually scores, and put it in front of checker.py and F's certified-faithful parse."""
import os,sys,json,pickle,collections,time
LAB='/home/user/integer_solver/solve_lab'; T=os.path.join(LAB,'agentT_work')
sys.path.insert(0,os.path.join(LAB,'agentM_work'))
sys.path.insert(0,os.path.join(LAB,'agentE_work'))
sys.set_int_max_str_digits(20_000_000)
os.chdir(os.path.join(LAB,'agentM_work'))
import ieng
D4=[642,28730,29854,31864]
t0=time.time()
r=ieng.tune(D4,want=True)
print('M engine on the witness subset %s'%D4)
print('   ok=%s  base_score=%s  score=%s  nknobs=%s  (%.0fs)'%(
      r.get('ok'),r.get('base_score'),r.get('score'),r.get('nknobs'),time.time()-t0),flush=True)
assert r.get('changes') is not None, 'engine returned no assignment'
v,aff=ieng.apply_delta(ieng.V_UNC,r['changes'],r['pin'])
print('   changed vars: %d ; downstream recomputed: %d'%(len(r['changes']),len(aff)),flush=True)
out=os.path.join(T,'t_meng_assign.json')
json.dump({'x_%d'%i:str(x) for i,x in enumerate(v) if x},open(out,'w'))
print('   wrote %s'%out,flush=True)
# --- independent check 1: checker.py's own loader/evaluator, in-process ---
sys.path.insert(0,LAB)
import checker as CK
codes,varsets=CK.load_equations()
fails=CK.evaluate_all(codes,v)
print('\nCHECKER (independent of M): satisfied %d/%d  (%d failing)'%(len(codes)-len(fails),len(codes),len(fails)),flush=True)
print('   failing: %s'%fails,flush=True)
DELIV=[12231,12270,12350,14584,18673,22044,29125]
print('   == the deliverable\'s exact 7 failures? %s'%(fails==DELIV),flush=True)
print('   M\'s engine reported %s ; checker says %d  -> AGREE: %s'%(
      r.get('score'),len(codes)-len(fails),r.get('score')==len(codes)-len(fails)),flush=True)
# --- independent check 2: F's certified-faithful atom parse ---
sys.path.insert(0,os.path.join(LAB,'agentF_work'))
from fwd import compile_node
d=pickle.load(open(os.path.join(LAB,'agentF_work','circ4.pkl'),'rb'))
atoms=d['atoms']; eqrows=d['eqrows']; names=list(atoms); idx={a:i for i,a in enumerate(names)}
a2e=collections.defaultdict(set)
for e,row in enumerate(eqrows):
    for k,a in row: a2e[idx[a]].add(e)
prog=compile('r[:]=['+','.join(compile_node(atoms[a]) for a in names)+']','<at>','exec')
rr=[0]*len(names); exec(prog,{'v':v,'r':rr,'__builtins__':{}})
nz=[i for i in range(len(names)) if rr[i]]
foot=set()
for i in nz: foot|=a2e[i]
print('\nF PARSE (independent of M): %d nonzero atoms'%len(nz),flush=True)
for i in nz: print('     %s'%names[i][:74],flush=True)
print('   equation footprint: %d ; == checker failing set? %s'%(len(foot),foot==set(fails)),flush=True)
