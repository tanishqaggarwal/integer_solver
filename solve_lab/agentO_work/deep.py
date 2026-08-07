"""Deeper closure: does raising maxr/maxv reach the pin repairs?  Reports blocking rows."""
import sys, json, os, time
ED='/home/user/integer_solver/solve_lab/agentE_work'
sys.path.insert(0,ED); os.chdir(ED)
sys.set_int_max_str_digits(20_000_000)
import channels as C, engine as E, fast, sparse
OD='/home/user/integer_solver/solve_lab/agentO_work'

def build(seed, maxr, maxv, verbose=True):
    v0=E.forward(seed); bad0=E.badatoms(v0)
    S=set(C.CLUSTERKN); pend=set(bad0); seenA=set(); cols={}; knobs=[]
    for rnd in range(maxr+1):
        new=set()
        for a in pend: new|=set(E.cone(a)[1])
        new-=S|set(seed)
        if not new: break
        t0=time.time()
        aff,c2=C.affine_cols(v0,bad0,sorted(new))
        cols.update(c2); knobs+=aff; S|=set(new)
        if verbose: print(f"   rnd{rnd}: cand={len(new)} affine={len(aff)} tot_knobs={len(knobs)} S={len(S)} {time.time()-t0:.0f}s",flush=True)
        t=set()
        for f in aff: t|=set(cols[f])
        seenA|=pend; pend=(t|set(bad0))-seenA
        if len(S)>maxv: break
    aff0,c0=C.affine_cols(v0,bad0,C.CLUSTERKN)
    for f in aff0:
        if f not in cols: knobs.append(f); cols[f]=c0[f]
    rows_at=set(bad0)
    for f in knobs: rows_at|=set(cols[f])
    rows_at=sorted(rows_at)
    rowmap={a:{} for a in rows_at}
    for f in knobs:
        for a,c in cols[f].items(): rowmap[a][f]=c
    return v0,bad0,knobs,cols,rows_at,rowmap

def solve(v0,bad0,knobs,cols,rows_at,rowmap,seed,verbose=True):
    sol,msg,_=sparse.solve_sparse([rowmap[a] for a in rows_at],[-bad0.get(a,0) for a in rows_at],
                                  names=rows_at,verbose=False,maxcore=400,maxcorebits=5_000_000)
    blocking=[]
    if sol is None:
        if verbose: print('   full system:',msg,flush=True)
        keep=[]
        for i,a in enumerate(rows_at):
            idx=keep+[i]
            s2,_,_=sparse.solve_sparse([rowmap[rows_at[j]] for j in idx],[-bad0.get(rows_at[j],0) for j in idx],
                                       verbose=False,maxcore=400,maxcorebits=5_000_000)
            if s2 is not None: keep=idx
            else: blocking.append(a)
        sol,_,_=sparse.solve_sparse([rowmap[rows_at[j]] for j in keep],[-bad0.get(rows_at[j],0) for j in keep],
                                    verbose=False,maxcore=400,maxcorebits=5_000_000)
        if verbose: print(f'   kept {len(keep)}/{len(rows_at)} blocking={blocking}',flush=True)
    ns=dict(seed)
    for f,d in (sol or {}).items():
        if d: ns[f]=v0[f]+d
    v=E.forward(ns); av=E.badatoms(v)
    return len(E.eqfails(av)),ns,sorted(av),v,blocking

if __name__=='__main__':
    bits=[int(x) for x in sys.argv[1].split(',')] if sys.argv[1]!='empty' else []
    maxr=int(sys.argv[2]); maxv=int(sys.argv[3])
    s=dict(C.base)
    for b in bits: s[b]=1
    print(f'=== bits={bits} maxr={maxr} maxv={maxv}',flush=True)
    t0=time.time()
    v0,bad0,knobs,cols,rows_at,rowmap=build(s,maxr,maxv)
    print(f'   bad0={sorted(bad0)} knobs={len(knobs)} rows={len(rows_at)} ({time.time()-t0:.0f}s)',flush=True)
    n,ns,av,v,blk=solve(v0,bad0,knobs,cols,rows_at,rowmap,s)
    print(f'   RESULT fails={n} score={39033-n} bad={av} blocking={blk} ({time.time()-t0:.0f}s)',flush=True)
