"""Decisive: what is the achievable set of (A mod p, B mod p)?
Relax atoms 22229 (=A), 22231 (=B) and the square 37887, making x_7068 and x_4432 free inputs.
All other checks must stay 0.  Then test which of the 9 cancellation ratios are reachable."""
import pickle, collections, heapq, sys, os, time
HERE=os.path.dirname(os.path.abspath(__file__)); S9=os.path.dirname(HERE)
sys.path.insert(0,S9); os.chdir(S9)
import harness as H, modsolve
P = 2**256-2**32-977
NV = 38748
g = pickle.load(open('gates.pkl','rb')); polys = pickle.load(open('polys.pkl','rb'))
topo = pickle.load(open('topo.pkl','rb'))['topo']; roots = pickle.load(open('roots.pkl','rb'))
definer = dict(g['definer']); atom_out = dict(g['atom_out'])
for a, vv in ((22229,7068),(22231,4432)):
    del atom_out[a]; del definer[vv]
RELAXED = {22229, 22231, 37887}
avars = [set(v for m in Pp for v in m) for Pp in polys]
var_atoms = collections.defaultdict(list)
for a, s in enumerate(avars):
    for u in s: var_atoms[u].append(a)
tidx = {v: i for i, v in enumerate(topo)}; BIG = len(topo)+10
checks = [a for a in range(len(polys)) if a not in atom_out and a not in RELAXED]
rp = {a: (roots[a] if a in roots else polys[a]) for a in checks}
freeinp = [x for x in range(NV) if x not in definer]
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
    old = v[t]; v[t] = 0; rest = evalpoly(Pp, v); v[t] = old
    if rest % c: return None
    return -rest//c
def ripple(v, seeds):
    pq = []
    for u, val in seeds.items():
        if v[u]==val: continue
        v[u]=val; heapq.heappush(pq,(tidx.get(u,BIG),u))
    while pq:
        _, u = heapq.heappop(pq)
        for a in var_atoms[u]:
            if a not in atom_out: continue
            c, t = atom_out[a]
            if t == u: continue
            nv = solve_for(a, t, v)
            if nv is None or nv == v[t]: continue
            v[t] = nv; heapq.heappush(pq,(tidx.get(t,BIG),t))

if __name__ == '__main__':
    v = H.load_assignment('../best/new_instance_partial_39022.json')
    base = {a: evalpoly(Pp, v) for a, Pp in rp.items()}
    print('nonzero checks at baseline (both chains relaxed):', {a:x for a,x in base.items() if x})
    t0=time.time(); J={}
    for i,f in enumerate(freeinp):
        vv=list(v); ripple(vv,{f:v[f]+1})
        col={}
        for a,Pp in rp.items():
            nv=evalpoly(Pp,vv)
            if nv!=base[a]: col[a]=nv-base[a]
        if col: J[f]=col
    print(f'jac {len(J)} cols, {time.time()-t0:.0f}s')
    pickle.dump({'J':J,'base':base}, open('chain2/jacAB.pkl','wb'))
    zero = {a:0 for a in base}
    for name,(al,be) in [('(1,0)',(1,0)), ('(0,1)',(0,1)), ('(1,1)',(1,1))]:
        extra = {'ROW_A': ({7068:1}, al), 'ROW_B': ({4432:1}, be)}
        sol,msg = modsolve.solve(J, zero, extra_rows=extra)
        print(f'  achievable delta (d_x7068, d_x4432) = {name}: {"YES" if sol is not None else "NO  ("+msg+")"}')
