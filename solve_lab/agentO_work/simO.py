"""Agent O simultaneous solve: E's LOG-16 idea but with the UNFILTERED knob set
   (bitfeas2-style: keep every cone-free knob, drop non-affine (knob,row) PAIRS by
   dropping the row, never the knob) + exact re-propagation iteration."""
import sys, json, os, time, collections
ED='/home/user/integer_solver/solve_lab/agentE_work'
sys.path.insert(0,ED)
if os.getcwd()!=ED: os.chdir(ED)
sys.set_int_max_str_digits(20_000_000)
import engine as E, fast, sparse, harness as H, channels as C
import re
OD='/home/user/integer_solver/solve_lab/agentO_work'

_isb={}
def isbool(f):
    r=_isb.get(f)
    if r is not None: return r
    r=False
    for i in H.occ[f]:
        t=re.sub(r'x_%d\b'%f,'X',H.atoms[i])
        if t in ('X - X * X','X * X - X','X * (X - 1)','2 * X * (1 - X)'): r=True; break
    _isb[f]=r; return r

def closure(v0, bad0, frozen, maxr, maxv, verbose=False):
    """Returns knobs (list), cols {f:{atom:delta}}, nonlin set of (f,atom), rounds."""
    S=[]; cols={}; nonlin=set(); processed=set(); rounds={}
    pending=set(bad0)
    for rnd in range(maxr+1):
        newS=set()
        for a in pending: newS|=set(E.cone(a)[1])
        newS-=set(S)|frozen
        newS=sorted(newS)
        if not newS: break
        t0=time.time()
        for f in newS:
            if isbool(f): continue          # booleans are 0/1 decisions, not knobs
            b1,_=fast.resid_delta(v0,bad0,{f:v0[f]+1})
            b2,_=fast.resid_delta(v0,bad0,{f:v0[f]+2})
            col={}
            for a in set(b1)|set(bad0):
                d=b1.get(a,0)-bad0.get(a,0)
                if d: col[a]=d
            for a in set(b2)|set(bad0)|set(col):
                if b2.get(a,0)-bad0.get(a,0)!=2*col.get(a,0): nonlin.add((f,a))
            cols[f]=col; S.append(f); rounds.setdefault(rnd,[]).append(f)
        aff=set()
        for f in newS:
            if f in cols: aff|=set(cols[f])
        processed|=pending; pending=(aff|set(bad0))-processed
        if verbose: print(f"    rnd{rnd}: cand={len(newS)} knobs+={len(rounds.get(rnd,[]))} tot={len(S)} {time.time()-t0:.0f}s",flush=True)
        if len(S)>maxv: break
    return S,cols,nonlin,rounds

def make_system(bad0,S,cols,nonlin,keep_dirty=False):
    Sset=set(S); atoms=set(bad0)
    for f in S: atoms|=set(cols[f])
    nl={a for f,a in nonlin if f in Sset}
    if keep_dirty: use=sorted(atoms)
    else: use=sorted(a for a in atoms if a not in nl)
    rowmap={a:{} for a in use}
    for f in S:
        for a,c in cols[f].items():
            if a in rowmap: rowmap[a][f]=c
    return use,rowmap,sorted(nl)

def maxsolvable(use,rowmap,bad0,prio=None,verbose=False,**kw):
    """Greedy maximal solvable subset, bad0 rows tried FIRST."""
    order=sorted(range(len(use)), key=lambda i:(0 if use[i] in bad0 else 1, i))
    keep=[]; blocked=[]
    for i in order:
        idx=keep+[i]
        s2,_,_=sparse.solve_sparse([rowmap[use[j]] for j in idx],[-bad0.get(use[j],0) for j in idx],
                                   verbose=False,**kw)
        if s2 is not None: keep=idx
        else: blocked.append(use[i])
    sol,msg,_=sparse.solve_sparse([rowmap[use[j]] for j in keep],[-bad0.get(use[j],0) for j in keep],
                                  verbose=False,**kw)
    return sol,keep,blocked

def solve_once(seed, frozen, maxr=5, maxv=4000, keep_dirty=False, verbose=False,
               maxcore=400, maxcorebits=5_000_000):
    v0=E.forward(seed); bad0=E.badatoms(v0)
    if not bad0: return 0,seed,[],v0,[],{}
    S,cols,nonlin,rounds=closure(v0,bad0,frozen,maxr,maxv,verbose=verbose)
    use,rowmap,nl=make_system(bad0,S,cols,nonlin,keep_dirty)
    info={'knobs':len(S),'rows':len(use),'dropped_nonlin':len(nl),'bad0':sorted(bad0)}
    if verbose: print('    system:',info,flush=True)
    kw=dict(maxcore=maxcore,maxcorebits=maxcorebits)
    sol,msg,_=sparse.solve_sparse([rowmap[a] for a in use],[-bad0.get(a,0) for a in use],
                                  names=use,verbose=False,**kw)
    blocked=[]
    if sol is None:
        info['full_msg']=msg
        if verbose: print('    full unsat:',msg,'-> greedy over',len(use),'rows',flush=True)
        sol,keep,blocked=maxsolvable(use,rowmap,bad0,**kw)
        info['kept']=len(keep)
    info['blocked']=blocked
    ns=dict(seed)
    for f,d in (sol or {}).items():
        if d: ns[f]=v0[f]+d
    v=E.forward(ns); av=E.badatoms(v)
    return len(E.eqfails(av)),ns,sorted(av),v,blocked,info

def iterate(seed, frozen, iters=5, **kw):
    s=dict(seed); best=None; hist=[]
    for it in range(iters):
        n,ns,av,v,blk,info=solve_once(s,frozen,**kw)
        hist.append((n,av,info))
        if best is None or n<best[0]: best=(n,dict(ns),av,v)
        if n==0: return 0,ns,av,v,hist
        if ns==s: break
        s=ns
    return best[0],best[1],best[2],best[3],hist
