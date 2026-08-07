"""SIMULTANEOUS composition: flip a bit whose pin system is known integrally solvable, and
   solve its pin rows TOGETHER with the cluster's affine knobs as ONE system (not sequentially)."""
import sys, json, math, pickle, time, re
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import engine as E, fast, sparse, harness as H
P=115792089237316195423570985008687907853269984665640564039457584007908834671663
CLUSTERKN=[6083,11436,14393,14853,22820,26489,31339,37012]
base={int(k):int(v) for k,v in json.load(open('triple8_seed.json')).items()}

def affine_cols(v0,bad0,cand):
    cols={}; aff=[]
    for f in cand:
        o=v0[f]
        b1,_=fast.resid_delta(v0,bad0,{f:o+1})
        b2,_=fast.resid_delta(v0,bad0,{f:o+2})
        b7,_=fast.resid_delta(v0,bad0,{f:o+7})
        col={}; ok=True
        for a in set(b1)|set(b2)|set(b7)|set(bad0):
            d1=b1.get(a,0)-bad0.get(a,0)
            if b2.get(a,0)-bad0.get(a,0)!=2*d1 or b7.get(a,0)-bad0.get(a,0)!=7*d1: ok=False; break
            if d1: col[a]=d1
        if ok: aff.append(f); cols[f]=col
    return aff,cols

def run(bit, maxr=4, maxv=2500, log=sys.stdout):
    s=dict(base); s[bit]=1
    v0=E.forward(s); bad0=E.badatoms(v0)
    print(f"[x_{bit}] flipped: bad={sorted(bad0)}",file=log,flush=True)
    # ONE closure over every bad atom (cluster rows AND the bit's pin rows), affine knobs only
    S=set(CLUSTERKN); pend=set(bad0); seenA=set()
    for rnd in range(maxr+1):
        new=set()
        for a in pend: new|=set(E.cone(a)[1])
        new-= S|{18956,1530,1603,bit}
        if not new: break
        aff,c2=affine_cols(v0,bad0,sorted(new))
        if rnd==0: cols=dict(c2); knobs=list(aff)
        else: cols.update(c2); knobs+=aff
        S|=set(new)
        touched=set()
        for f in aff: touched|=set(cols[f])
        seenA|=pend; pend=(touched|set(bad0))-seenA
        if len(S)>maxv: break
    aff0,c0=affine_cols(v0,bad0,CLUSTERKN)
    for f in aff0:
        if f not in cols: knobs.append(f); cols[f]=c0[f]
    rows_at=set(bad0)
    for f in knobs: rows_at|=set(cols[f])
    rows_at=sorted(rows_at)
    print(f"   ONE system: {len(knobs)} affine knobs x {len(rows_at)} rows",file=log,flush=True)
    rowmap={a:{} for a in rows_at}
    for f in knobs:
        for a,c in cols[f].items(): rowmap[a][f]=c
    sol,msg,_=sparse.solve_sparse([rowmap[a] for a in rows_at],[-bad0.get(a,0) for a in rows_at],
                                  names=rows_at,verbose=False,maxcore=500,maxcorebits=5_000_000)
    print(f"   -> {msg[:90]}",file=log,flush=True)
    if sol is None:
        keep=[]
        for i,a in enumerate(rows_at):
            idx=keep+[i]
            s2,_,_=sparse.solve_sparse([rowmap[rows_at[j]] for j in idx],[-bad0.get(rows_at[j],0) for j in idx],
                                       verbose=False,maxcore=500,maxcorebits=5_000_000)
            if s2 is not None: keep=idx
        drop=[rows_at[i] for i in range(len(rows_at)) if i not in keep]
        print(f"   max solvable {len(keep)}/{len(rows_at)}  BLOCKING {drop}",file=log,flush=True)
        sol,_,_=sparse.solve_sparse([rowmap[rows_at[j]] for j in keep],[-bad0.get(rows_at[j],0) for j in keep],
                                    verbose=False,maxcore=500,maxcorebits=5_000_000)
    if sol is None: return None
    ns=dict(s)
    for f,d in sol.items():
        if d: ns[f]=v0[f]+d
    v=E.forward(ns); av=E.badatoms(v); ff=E.eqfails(av)
    print(f"   EXACT fails={len(ff)} score={39033-len(ff)} bad={sorted(av)[:12]}",file=log,flush=True)
    return len(ff),ns,sorted(av),v

if __name__=='__main__':
    bits=[int(x) for x in sys.argv[1:]]
    best=None
    for b in bits:
        t0=time.time()
        try: r=run(b)
        except Exception as e:
            print(f"[x_{b}] ERR {type(e).__name__} {e}",flush=True); continue
        print(f"[x_{b}] done {time.time()-t0:.0f}s",flush=True)
        if r and (best is None or r[0]<best[0]):
            best=r
            json.dump({f"x_{i}":int(r[3][i]) for i in range(E.NV) if r[3][i]!=0}, open('simul_%d.json'%(39033-r[0]),'w'))
            json.dump({str(k):str(int(x)) for k,x in r[1].items()}, open('simul_seed.json','w'))
    print("BEST",best[0] if best else None, best[2] if best else None,flush=True)
