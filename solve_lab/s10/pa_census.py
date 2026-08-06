"""WHOLE-INSTANCE break census.

For every variable t we perturb it away from the delivered witness while BLOCKING its
definer atom (so that atom stays broken) and let the gate ripple repair everything else.
The result is the exact atom support S reachable by "breaking atom definer[t]", measured
by construction, together with |E(S)| and the true failing count.

Usage:  python3 pa_census.py <lo> <hi> <outfile>
"""
import os, sys, collections, json, time, random
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L

BROKEN_GATES={22229,22230,35758,35761,35762}
v0=L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
base_av=L.all_atom_values(v0)
base_nz=set(a for a in range(L.NA) if base_av[a])
base_fail=set(L.failing_eqs(base_av))
print('base nz',sorted(base_nz),'fail',len(base_fail),flush=True)
random.seed(11)
DELTA=random.randrange(1,2**80)

def probe(t, delta):
    w=list(v0)
    blk=set(BROKEN_GATES)
    d=L.definer.get(t)
    if d is not None: blk.add(d)
    try:
        changed,steps=L.ripple(w,{t:v0[t]+delta},maxsteps=60000,block=blk)
    except Exception as e:
        return None
    cand=set()
    for u in changed: cand.update(L.var_atoms[u])
    nz=set(base_nz); diff=set()
    for a in cand:
        nv=L.evalpoly(L.polys[a],w)
        if nv!=base_av[a]: diff.add(a)
        if nv: nz.add(a)
        else: nz.discard(a)
    # failing equations: recheck all equations touched by diff atoms, plus base fails
    eqs=set(base_fail)
    for a in diff: eqs.update(L.atom2eq.get(a,()))
    av={}
    fail=0; faillist=[]
    for e in eqs:
        m,sq,co=L.eq_atoms[e]
        s=0
        for a,c in co.items():
            if a not in av: av[a]=L.evalpoly(L.polys[a],w)
            s+=c*av[a]
        if s: fail+=1; faillist.append(e)
    E=set()
    for a in nz: E.update(L.atom2eq.get(a,()))
    return dict(nz=sorted(nz),nE=len(E),fail=fail,faillist=sorted(faillist)[:20],steps=steps,nch=len(changed))

if __name__=='__main__':
    lo,hi,out=int(sys.argv[1]),int(sys.argv[2]),sys.argv[3]
    res={}
    t0=time.time()
    for t in range(lo,min(hi,L.NVARS)):
        r=probe(t,DELTA)
        if r is None: continue
        r['a']=L.definer.get(t)
        res[t]=r
        if (t-lo)%200==0:
            print(f'{t} fail={r["fail"]} nz={len(r["nz"])} nE={r["nE"]} ch={r["nch"]} {time.time()-t0:.0f}s',flush=True)
    json.dump(res,open(out,'w'))
    best=sorted(res.items(),key=lambda kv:kv[1]['fail'])[:25]
    print('BEST in range:')
    for t,r in best:
        print(f'  x_{t} (definer a{r["a"]}): fail={r["fail"]} |S|={len(r["nz"])} |E|={r["nE"]} S={r["nz"][:12]}')
    print(f'{time.time()-t0:.0f}s done')
