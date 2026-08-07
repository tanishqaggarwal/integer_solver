"""Channel model: measure channels at a configuration, cross-check the partition, then
   enumerate channel-set x representative EXACTLY with the simultaneous solve of LOG 16."""
import sys, json, re, collections, itertools, time, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import engine as E, fast, sparse, harness as H
P=115792089237316195423570985008687907853269984665640564039457584007908834671663
ROWS=[7389,10187,20212,20215,28647]
CLUSTERKN=[6083,11436,14393,14853,22820,26489,31339,37012]
base={int(k):int(v) for k,v in json.load(open('/home/user/integer_solver/solve_lab/agentE_work/triple8_seed.json')).items()}
def isb(f):
    for i in H.occ[f]:
        t=re.sub(r'x_%d\b'%f,'X',H.atoms[i])
        if t in ('X - X * X','X * X - X','X * (X - 1)','2 * X * (1 - X)'): return True
    return False
CAND=sorted(set().union(*[set(E.cone(a)[1]) for a in ROWS]))
BOOLS=[f for f in CAND if isb(f)]

def channels(seed):
    v0=E.forward(seed); bad0=E.badatoms(v0)
    co={}
    for f in BOOLS:
        if v0[f]==1: continue
        b1,_=fast.resid_delta(v0,bad0,{f:1})
        d={a:(b1.get(a,0)-bad0.get(a,0)) for a in ROWS}
        c=((d[20212]+d[28647])%P,(d[20215]+d[10187])%P)
        if any(c): co[f]=c
    cls=collections.defaultdict(list)
    for f,c in co.items(): cls[c].append(f)
    return v0,bad0,cls

def affine_cols(v0,bad0,cand):
    cols={}; aff=[]
    for f in cand:
        o=v0[f]
        b1,_=fast.resid_delta(v0,bad0,{f:o+1}); b2,_=fast.resid_delta(v0,bad0,{f:o+2}); b7,_=fast.resid_delta(v0,bad0,{f:o+7})
        col={}; ok=True
        for a in set(b1)|set(b2)|set(b7)|set(bad0):
            d1=b1.get(a,0)-bad0.get(a,0)
            if b2.get(a,0)-bad0.get(a,0)!=2*d1 or b7.get(a,0)-bad0.get(a,0)!=7*d1: ok=False; break
            if d1: col[a]=d1
        if ok: aff.append(f); cols[f]=col
    return aff,cols

def simsolve(seed, maxr=3, maxv=2000):
    """LOG 16 simultaneous solve at this seed."""
    v0=E.forward(seed); bad0=E.badatoms(v0)
    if not bad0: return 0,seed,[],v0
    S=set(CLUSTERKN); pend=set(bad0); seenA=set(); cols={}; knobs=[]
    for rnd in range(maxr+1):
        new=set()
        for a in pend: new|=set(E.cone(a)[1])
        new-=S|set(seed)
        if not new: break
        aff,c2=affine_cols(v0,bad0,sorted(new))
        cols.update(c2); knobs+=aff; S|=set(new)
        t=set()
        for f in aff: t|=set(cols[f])
        seenA|=pend; pend=(t|set(bad0))-seenA
        if len(S)>maxv: break
    aff0,c0=affine_cols(v0,bad0,CLUSTERKN)
    for f in aff0:
        if f not in cols: knobs.append(f); cols[f]=c0[f]
    rows_at=set(bad0)
    for f in knobs: rows_at|=set(cols[f])
    rows_at=sorted(rows_at)
    rowmap={a:{} for a in rows_at}
    for f in knobs:
        for a,c in cols[f].items(): rowmap[a][f]=c
    sol,msg,_=sparse.solve_sparse([rowmap[a] for a in rows_at],[-bad0.get(a,0) for a in rows_at],
                                  names=rows_at,verbose=False,maxcore=400,maxcorebits=5_000_000)
    if sol is None:
        keep=[]
        for i,a in enumerate(rows_at):
            idx=keep+[i]
            s2,_,_=sparse.solve_sparse([rowmap[rows_at[j]] for j in idx],[-bad0.get(rows_at[j],0) for j in idx],
                                       verbose=False,maxcore=400,maxcorebits=5_000_000)
            if s2 is not None: keep=idx
        sol,_,_=sparse.solve_sparse([rowmap[rows_at[j]] for j in keep],[-bad0.get(rows_at[j],0) for j in keep],
                                    verbose=False,maxcore=400,maxcorebits=5_000_000)
    if sol is None: return None
    ns=dict(seed)
    for f,d in sol.items():
        if d: ns[f]=v0[f]+d
    v=E.forward(ns); av=E.badatoms(v)
    return len(E.eqfails(av)),ns,sorted(av),v
