"""Iterative goal-directed repair: forward-eval, then solve each bad check for a free input."""
import sys, os, json, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw

P=L.P
FREE=[u for u in range(L.NVARS) if L.definer.get(u) is None]
FREESET=set(FREE)
NAT={u:len(L.var_atoms[u]) for u in range(L.NVARS)}

def free_cands(a):
    """free vars in atom a, occurring linearly, fewest-atoms first"""
    out=[]
    for u in L.avars[a]:
        if u not in FREESET: continue
        ok=True
        for m in L.polys[a]:
            if m.count(u)>1: ok=False; break
        if ok: out.append(u)
    out.sort(key=lambda u: (NAT[u], u))
    return out

def run(sel, rounds=60, verbose=True):
    v=[0]*L.NVARS
    for k,x in sel.items(): v[k]=x
    locked=set(sel)
    hist=[]
    for it in range(rounds):
        fw.forward(v)
        bad=fw.bad_checks(v)
        av=L.all_atom_values(v); f=L.failing_eqs(av)
        hist.append((len(bad),len(f)))
        if verbose: print(f"  it{it}: bad_checks={len(bad)} failing_eqs={len(f)} score={L.NEQ-len(f)}")
        if not bad: return v, bad, f
        prog=False
        for a in bad:
            if L.polys[a] and fw.evalpoly(L.polys[a],v)==0: continue
            for t in free_cands(a):
                if t in locked: continue
                x=fw.solve_lin(a,t,v)
                if x is not None and x!=v[t]:
                    v[t]=x; prog=True; break
        if not prog:
            if verbose: print("  no progress")
            break
    fw.forward(v)
    bad=fw.bad_checks(v); av=L.all_atom_values(v); f=L.failing_eqs(av)
    return v, bad, f

if __name__=='__main__':
    sel={542:1, 438:1}
    t0=time.time()
    v,bad,f=run(sel)
    print(f"FINAL bad_checks={len(bad)} failing={len(f)} score={L.NEQ-len(f)}  ({time.time()-t0:.0f}s)")
    print("bad:", bad[:30])
    json.dump({str(i):v[i] for i in range(L.NVARS)}, open('rep_state.json','w'))
