"""Deep repair: find free vars in a check's CONE with an exact linear effect on it."""
import sys, os, json, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw

def cone(a):
    seen=set(); st=[u for u in L.avars[a]]
    while st:
        u=st.pop()
        if u in seen: continue
        seen.add(u)
        d=L.definer.get(u)
        if d is None: continue
        for w in L.avars[d]:
            if w!=u: st.append(w)
    return seen

def cone_free(a):
    return sorted(u for u in cone(a) if L.definer.get(u) is None)

def sub_order(cs):
    """restricted forward order: SCC comps intersecting the cone"""
    return [c for c in fw.ORDER if any(u in cs for u in c)]

def local_forward(v, order):
    for comp in order:
        if len(comp)==1:
            u=comp[0]
            x=fw.solve_lin(L.definer[u], u, v)
            if x is not None: v[u]=x
        else:
            for _ in range(40):
                ch=False
                for u in comp:
                    x=fw.solve_lin(L.definer[u], u, v)
                    if x is not None and x!=v[u]: v[u]=x; ch=True
                if not ch: break

def handles(v, a, locked=(), limit=None):
    """return [(var, delta_per_unit)] for free vars with exact linear effect on atom a"""
    cs=cone(a); order=sub_order(cs)
    frees=[u for u in sorted(cs) if L.definer.get(u) is None and u not in locked]
    if limit: frees=frees[:limit]
    base=fw.evalpoly(L.polys[a], v)
    out=[]
    snap={u:v[u] for u in cs}
    for t in frees:
        o=v[t]
        v[t]=o+1; local_forward(v, order); d1=fw.evalpoly(L.polys[a],v)-base
        for u,val in snap.items(): v[u]=val
        v[t]=o
        if d1==0: continue
        v[t]=o+2; local_forward(v, order); d2=fw.evalpoly(L.polys[a],v)-base
        for u,val in snap.items(): v[u]=val
        v[t]=o
        if d2==2*d1: out.append((t,d1))
    out.sort(key=lambda kv:(len(L.var_atoms[kv[0]]), kv[0]))
    return out, base

if __name__=='__main__':
    v=[int(x) for _,x in sorted(((int(k),val) for k,val in json.load(open('solve.json')).items()))]
    fw.forward(v)
    for a in [26719,26721,26723,26733,28438,32342,36185]:
        t0=time.time()
        h,base=handles(v,a)
        ok=[(t,d) for t,d in h if base % d==0]
        print(f"a{a}: base={str(base)[:28]}... conefree={len(cone_free(a))} linear_handles={len(h)} exact={len(ok)}  ({time.time()-t0:.0f}s)")
        print("    top:", [(t,d) for t,d in h[:6]])
