import pickle, subprocess, sys, os, collections
W='/home/user/integer_solver/solve_lab/agentB_work/'
M=pickle.load(open(W+'model5.pkl','rb')); facs=M['facs']
fvars=[]
for p in facs:
    vs=set()
    for m in p: vs.update(m)
    fvars.append(vs)
import importlib.util
spec=importlib.util.spec_from_file_location('bo7', W+'borient7.py'); bo7=importlib.util.module_from_spec(spec)
import sys as _s; _s.argv=['x']; spec.loader.exec_module(bo7) if False else None
cands=[]
for p in facs:
    sq=set()
    for m in p:
        if len(m)==2 and m[0]==m[1]: sq.add(m[0])
    allv=set()
    for m in p: allv.update(m)
    ok=allv-sq
    unit=set(m[0] for m,c in p.items() if len(m)==1 and abs(c)==1 and m[0] in ok)
    quad=set()
    for m in p:
        if len(m)>1: quad.update(m)
    unit |= set(v for v in ok&quad if v not in [m[0] for m in p if len(m)==1])
    cands.append(unit if unit else ok)
forbid=set()
best=None
for rnd in range(12):
    pickle.dump(sorted(forbid), open(W+'forbid.pkl','wb'))
    r=subprocess.run([sys.executable, W+'borient7.py', W+'excl_base.pkl','orient7.pkl'],capture_output=True,text=True,cwd=W)
    line=[l for l in r.stdout.split('\n') if l.startswith('DEFINED')]
    O=pickle.load(open(W+'orient7.pkl','rb'))
    free=set(O['free']); ass=O['assertions']; excl=set(O['excl'])
    new=set()
    for f in ass:
        if f in excl: continue
        U=[v for v in fvars[f] if v in free and v in cands[f]]
        if len(U)==1: new.add(U[0])
    print('round %d: %s  convertible-blocked vars=%d  forbid=%d'%(rnd,line[0] if line else r.stdout[-200:],len(new),len(forbid)),flush=True)
    if not (new-forbid): break
    forbid |= new
    os.replace(W+'orient7.pkl', W+'orient7_r%d.pkl'%rnd)
    os.replace(W+'orient7_r%d.pkl'%rnd, W+'orient7.pkl')
