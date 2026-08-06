"""EQUATION-SPACE solve: atoms need NOT vanish individually, only their per-equation
combinations must.  Strictly weaker than the atom-space system used everywhere so far."""
import pickle, collections, heapq, sys, os, time
HERE=os.path.dirname(os.path.abspath(__file__)); S9=os.path.dirname(HERE)
sys.path.insert(0,S9); os.chdir(S9)
import harness as H, modsolve
P = 2**256-2**32-977
NV = 38748
g = pickle.load(open('gates.pkl','rb')); polys = pickle.load(open('polys.pkl','rb'))
topo = pickle.load(open('topo.pkl','rb'))['topo']; roots = pickle.load(open('roots.pkl','rb'))
d = pickle.load(open('atoms.pkl','rb')); eq_terms = d['eq_terms']
definer = dict(g['definer']); atom_out = dict(g['atom_out'])
for a, vv in ((22229,7068),(22231,4432)):
    del atom_out[a]; del definer[vv]
avars = [set(v for m in Pp for v in m) for Pp in polys]
var_atoms = collections.defaultdict(list)
for a, s in enumerate(avars):
    for u in s: var_atoms[u].append(a)
tidx = {v: i for i, v in enumerate(topo)}; BIG = len(topo)+10
freeinp = [x for x in range(NV) if x not in definer]
NONGATE = [a for a in range(len(polys)) if a not in atom_out]     # everything not forced by a gate
SQ = set(roots)                                                    # square atoms: track their ROOT
def evalpoly(Pp, v):
    s = 0
    for m, c in Pp.items():
        t = c
        for u in m: t *= v[u]
        s += t
    return s
def solve_for(a, t, v):
    Pp = polys[a]; c = 0
    for m, cc in Pp.items():
        if len(m)==1 and m[0]==t: c += cc
        elif t in m: return None
    if c == 0: return None
    old=v[t]; v[t]=0; rest=evalpoly(Pp,v); v[t]=old
    if rest % c: return None
    return -rest//c
def ripple(v, seeds):
    pq=[]
    for u,val in seeds.items():
        if v[u]==val: continue
        v[u]=val; heapq.heappush(pq,(tidx.get(u,BIG),u))
    while pq:
        _,u = heapq.heappop(pq)
        for a in var_atoms[u]:
            if a not in atom_out: continue
            c,t = atom_out[a]
            if t==u: continue
            nv = solve_for(a,t,v)
            if nv is None or nv==v[t]: continue
            v[t]=nv; heapq.heappush(pq,(tidx.get(t,BIG),t))

if __name__ == '__main__':
    v = H.load_assignment('../best/new_instance_partial_39022.json')
    # atom values (squares tracked by ROOT so the row is linear)
    def atomvals(vv): return {a: evalpoly(roots[a] if a in SQ else polys[a], vv) for a in NONGATE}
    base = atomvals(v)
    nzb = {a:x for a,x in base.items() if x}
    print('nonzero non-gate atoms at baseline:', nzb.keys())
    t0=time.time(); J={}
    for f in freeinp:
        vv=list(v); ripple(vv,{f:v[f]+1})
        col={}
        for a in NONGATE:
            nv=evalpoly(roots[a] if a in SQ else polys[a], vv)
            if nv!=base[a]: col[a]=nv-base[a]
        if col: J[f]=col
    print(f'atom jac: {len(J)} cols nnz={sum(len(c) for c in J.values())} {time.time()-t0:.0f}s')
    # ---- equation-space rows: S_i = sum_a c_ia * y_a   (squares contribute root^2 -> 0 to 1st order)
    eqrow = collections.defaultdict(dict); eqrhs = {}
    atom_eqs = collections.defaultdict(list)
    for i,(m,sq,tl) in enumerate(eq_terms):
        for c,a in tl: atom_eqs[a].append((i,c))
    touched = set()
    for f,col in J.items():
        for a,dv in col.items():
            if a in SQ: continue                     # square handled by its own root row
            for i,c in atom_eqs.get(a,[]):
                eqrow[i][f] = (eqrow[i].get(f,0) + c*dv) % P
                touched.add(i)
    for i in touched:
        m,sq,tl = eq_terms[i]
        eqrhs[i] = (-sum(c*base.get(a,0) for c,a in tl)) % P
    print(f'equation-space rows: {len(eqrow)}')
    # ---- square-root rows (keep roots at 0: conservative)
    extra = {}
    for a in SQ:
        row = {f: J[f][a] % P for f in J if a in J[f]}
        if row: extra[('SQ',a)] = (row, (-base[a]) % P)
    print(f'square-root rows: {len(extra)}')
    rowsJ = {}       # convert eqrow (row-major) into modsolve's column-major J form
    colmaj = collections.defaultdict(dict)
    for i,row in eqrow.items():
        for f,val in row.items():
            if val: colmaj[f][('EQ',i)] = val
    basemap = {('EQ',i): (-eqrhs[i]) % P for i in eqrhs}
    for key,(row,r) in extra.items():
        for f,val in row.items(): colmaj[f][key] = val
        basemap[key] = (-r) % P
    print(f'combined system: {len(colmaj)} cols, {len(basemap)} rows')
    sol,msg = modsolve.solve(dict(colmaj), basemap)
    print('EQUATION-SPACE SOLVE:', msg)
    if sol is not None:
        nz={c:x for c,x in sol.items() if x}
        print(f'{len(nz)} nonzero deltas mod p'); pickle.dump(nz, open('kernel/eqsol.pkl','wb'))
