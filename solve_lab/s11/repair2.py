"""Repair loop with lookahead, incremental atom re-evaluation, and a check-aware cost.

Breaking a GATE atom is free (the ripple repairs it downstream); breaking a CHECK atom is real
damage.  So score each candidate move by the number of broken atoms it leaves after a ripple,
and take the best.  Ties broken by preferring moves that touch fewer checks.
"""
import sys, os, json, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw, atomval
P=L.P
HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
IS_CHECK=[L.atom_out.get(a) is None for a in range(L.NA)]

def full_atoms(v):
    return [L.evalpoly(L.polys[a],v) for a in range(L.NA)]

def lin_solve(a,t,v):
    c=0
    for m,cc in L.polys[a].items():
        k=m.count(t)
        if k==0: continue
        if k>1: return None
        term=cc
        for u in m:
            if u!=t: term*=v[u]
        c+=term
    if c==0: return None
    old=v[t]; v[t]=0
    rest=L.evalpoly(L.polys[a],v)
    v[t]=old
    if rest%c: return None
    return -rest//c

def apply_move(v, AV, broken, t, x):
    """set v[t]=x, ripple, update AV/broken incrementally.  returns undo info"""
    snap={t:v[t]}
    ch,_=L.ripple(v, {t:x})
    snap.update({u:None for u in ch if u not in snap})   # placeholder
    return ch

def run(src, out, rounds=200, verbose=True):
    v=load_raw(src)
    AV=full_atoms(v)
    broken=set(a for a in range(L.NA) if AV[a]!=0)
    print(f"start broken atoms={len(broken)} {sorted(broken)}")
    hist=[]
    for it in range(rounds):
        if not broken:
            break
        cands=[]
        for a in sorted(broken):
            for t in sorted(L.avars[a]):
                d=L.definer.get(t)
                if d is not None and d!=a:
                    continue          # setting a gate output directly just breaks its own gate
                x=lin_solve(a,t,v)
                if x is None or x==v[t]: continue
                cands.append((a,t,x))
        best=None
        for a,t,x in cands:
            v2=list(v)
            L.ripple(v2,{t:x})
            touched=set()
            for u in range(L.NVARS):
                pass
            # incremental: recompute atoms of changed vars + currently broken
            AV2=dict()
            diff=[u for u in range(L.NVARS) if v2[u]!=v[u]]
            cand_atoms=set(broken)
            for u in diff: cand_atoms|=set(L.var_atoms[u])
            nb=set(a2 for a2 in broken if a2 not in cand_atoms)
            for a2 in cand_atoms:
                if L.evalpoly(L.polys[a2],v2)!=0: nb.add(a2)
            key=(len(nb), sum(1 for a2 in nb if IS_CHECK[a2]))
            if best is None or key<best[0]:
                best=(key,a,t,x,v2,nb)
        if best is None:
            print(f"  it{it}: STUCK with {len(broken)} broken: {sorted(broken)}"); break
        key,a,t,x,v2,nb = best
        if key[0]>=len(broken) and (a,t) in hist[-4:]:
            print(f"  it{it}: no improving move (best leaves {key[0]})"); break
        hist.append((a,t))
        v=v2; broken=nb
        if verbose:
            print(f"  it{it}: a{a} <- x{t}   broken={len(broken)} ({sum(1 for z in broken if IS_CHECK[z])} checks) {sorted(broken)[:10]}", flush=True)
        if key[0]>=len(hist) and it>6 and len(set(hist[-6:]))<=2:
            print("  cycling"); break
    AV=full_atoms(v)
    F=L.failing_eqs(AV)
    print(f"RESULT broken atoms={len([a for a in range(L.NA) if AV[a]!=0])} failing eqs={len(F)} score={L.NEQ-len(F)}")
    json.dump({('x_%d'%i):v[i] for i in range(L.NVARS)}, open(out,'w'))
    print("saved",out)
    return v

if __name__=='__main__':
    src=sys.argv[1] if len(sys.argv)>1 else os.path.join(HERE,'data','fix2_round.json')
    run(src, os.path.join(HERE,'data','repair2_out.json'))
