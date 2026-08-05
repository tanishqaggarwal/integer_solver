#!/usr/bin/env python3
"""Value-driven minimal re-orientation repair (custom, no SAT/SMT).

Keep best's gate orientation everywhere. Force x_9770<-27973 and x_3183<-27978
(their true primitive gates). Recompute forward; any primitive gate broken by a
value change must grab a new output variable (augmenting path toward a free input
or a redundant combination atom). Iterate until no primitive gate is violated."""
import json, sys, time
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
sys.setrecursionlimit(2_000_000)
NV_PRIM = 4

def main():
    t0 = time.time()
    atoms = load_atoms(); avars = [atom_vars(p) for p in atoms]
    prim = [len(avars[a]) <= NV_PRIM for a in range(len(atoms))]
    # linear-output candidates per atom
    outv = {}
    for a in range(len(atoms)):
        bad = set(); lin = set()
        for m in atoms[a]:
            if len(m) == 1: lin.add(m[0])
            else: bad.update(m)
        outv[a] = lin - bad
    cand_atoms = defaultdict(list)   # var -> primitive atoms that can output it
    for a in range(len(atoms)):
        if prim[a]:
            for v in outv[a]: cand_atoms[v].append(a)

    best = json.load(open('best/best_partial_39013.json'))
    val = [0]*NVARS
    for k, x in best.items(): val[int(k[2:])] = x

    prov = json.load(open('eval_order.json'))['prov']
    df = [None]*NVARS      # var -> defining atom
    owner = [None]*len(atoms)
    for v, p in enumerate(prov):
        if p and p[0] >= 0:
            df[v] = p[0]; owner[p[0]] = v

    def force(atom, var):
        """Make `atom` define `var` via augmenting; returns True on success."""
        # standard alternating DFS; endpoints: free var (df None) or combo owner
        def dfs(a, seen):
            for v in outv[a]:
                if v in seen: continue
                seen.add(v)
                cur = df[v]
                if cur is None:                      # v is free -> take it
                    df[v] = a; owner[a] = v; return True
                if cur == a: continue
                if not prim[cur]:                    # displace redundant combo
                    owner[cur] = None; df[v] = a; owner[a] = v; return True
                if dfs(cur, seen):                   # re-home the primitive cur
                    df[v] = a; owner[a] = v; return True
            return False
        # release atom's current var first
        if owner[atom] is not None:
            df[owner[atom]] = None; owner[atom] = None
        return dfs(atom, set([var])) or _take(atom, var)
    def _take(atom, var):
        if df[var] is None or not prim[df[var]]:
            if df[var] is not None: owner[df[var]] = None
            df[var] = atom; owner[atom] = var; return True
        return False

    # force the two roots
    for v, a in [(9770, 27973), (3183, 27978)]:
        if owner[a] is not None and owner[a] != v:
            df[owner[a]] = None; owner[a] = None
        if df[v] is not None: owner[df[v]] = None
        df[v] = a; owner[a] = v

    order = json.load(open('eval_order.json'))['order']
    def gate_terms(v):
        a = df[v]; coef = 0; terms = []
        for m, c in atoms[a].items():
            if m == (v,): coef += c
            else: terms.append((c, m))
        return coef, terms
    def forward():
        seq = [v for v in order if df[v] is not None] + [v for v in range(NVARS) if df[v] is not None and v not in set(order)]
        for _ in range(40):
            ch = 0
            for v in seq:
                coef, terms = gate_terms(v); rs = 0
                for c, mv in terms:
                    t = c
                    for x in mv: t *= val[x]
                    rs += t
                if coef and (-rs) % coef == 0:
                    nv = (-rs)//coef
                    if nv != val[v]: val[v] = nv; ch += 1
            if ch == 0: break
    def viol_prims():
        out = []
        for a in range(len(atoms)):
            if not prim[a]: continue
            s = 0
            for m, c in atoms[a].items():
                t = c
                for x in m: t *= val[x]
                s += t
            if s: out.append(a)
        return out

    for it in range(30):
        forward()
        vp = viol_prims()
        if not vp:
            print(f"iter {it}: no violated primitives ({time.time()-t0:.0f}s)", flush=True); break
        fixed = 0
        for a in vp:
            # try to give atom a an output it can satisfy
            if owner[a] is not None: continue
            got = False
            for v in outv[a]:
                if force(a, v): got = True; break
            if got: fixed += 1
        print(f"iter {it}: violated prims {len(vp)}, re-homed {fixed} ({time.time()-t0:.0f}s)", flush=True)
        if fixed == 0: break

    forward()
    viol = []
    for a in range(len(atoms)):
        s = 0
        for m, c in atoms[a].items():
            t = c
            for x in m: t *= val[x]
            s += t
        if s: viol.append(a)
    print(f"FINAL violated atoms: {len(viol)}  {viol[:40]}", flush=True)
    json.dump({f"x_{i}": val[i] for i in range(NVARS)}, open('cand_augrepair.json', 'w'))
    print("wrote cand_augrepair.json", flush=True)

if __name__ == '__main__':
    main()
