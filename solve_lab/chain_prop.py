#!/usr/bin/env python3
"""Propagate the x_9770 / x_3183 correction through the primitive-gate chain.

Starting from the roots' true values, repeatedly: for each primitive gate that a
changed variable feeds, if it is now violated, solve it for its other linear
variable and record that as the new target. This follows the identity/sum chain
(x_9770 -> x_18274 -> x_6773 -> x_26517 -> ...) outward until it terminates at
combination atoms (redundant checks). Then apply and run the full checker."""
import json, time
from collections import deque, defaultdict
from propagate import load_atoms, atom_vars, NVARS
NV_PRIM = 4

def main():
    t0 = time.time()
    atoms = load_atoms(); avars = [atom_vars(p) for p in atoms]
    prim = [len(avars[a]) <= NV_PRIM for a in range(len(atoms))]
    var_atoms = defaultdict(list)
    for a in range(len(atoms)):
        if prim[a]:
            for v in avars[a]: var_atoms[v].append(a)
    # product-output map: v -> list of operand-sets (gates where v = product of operands)
    prod_out = defaultdict(list)
    for a in range(len(atoms)):
        if not prim[a]: continue
        lin = set(); bad = set()
        for m in atoms[a]:
            if len(m) == 1: lin.add(m[0])
            else: bad.update(m)
        for v in (lin - bad):
            ops = set()
            for m in atoms[a]:
                if len(m) >= 2: ops.update(m)
            if ops: prod_out[v].append(ops)

    best = json.load(open('best/best_partial_39019.json'))
    val = [0]*NVARS
    for k, x in best.items(): val[int(k[2:])] = x

    # targets: vars we've re-derived (already includes 9770,3183 in 39019)
    targets = set([9770, 3183])
    dq = deque([9770, 3183])
    stuck = []
    def solve_for(a, u):
        coef = 0; rest = 0
        for m, c in atoms[a].items():
            if m == (u,): coef += c
            else:
                t = c
                for x in m: t *= val[x]
                rest += t
        if coef == 0 or (-rest) % coef != 0: return None
        return (-rest) // coef
    while dq:
        v = dq.popleft()
        for a in var_atoms[v]:
            # is a violated now?
            s = 0
            for m, c in atoms[a].items():
                t = c
                for x in m: t *= val[x]
                s += t
            if s == 0: continue
            # find a linear var not yet a target, solvable
            lin = set(); bad = set()
            for m in atoms[a]:
                if len(m) == 1: lin.add(m[0])
                else: bad.update(m)
            cands = [u for u in (lin - bad) if u not in targets]
            # exclude vars pinned by a product gate whose operands are all unchanged
            def pinned(u):
                for ops in prod_out.get(u, []):
                    if not (ops & targets):  # this product fixes u to best
                        return True
                return False
            ranked = sorted(cands, key=lambda u: (pinned(u), len(var_atoms[u])))
            cands = [u for u in ranked if not pinned(u)] or ranked
            done = False
            for u in cands:
                nv = solve_for(a, u)
                if nv is not None:
                    val[u] = nv; targets.add(u); dq.append(u); done = True; break
            if not done:
                stuck.append(a)
    print(f"chain propagated: {len(targets)} vars re-derived, {len(set(stuck))} stuck primitive atoms ({time.time()-t0:.0f}s)", flush=True)
    print(f"  re-derived (first 40): {sorted(targets)[:40]}", flush=True)

    # violated atoms now
    viol = []
    for a in range(len(atoms)):
        s = 0
        for m, c in atoms[a].items():
            t = c
            for x in m: t *= val[x]
            s += t
        if s: viol.append(a)
    vp = [a for a in viol if prim[a]]
    print(f"violated atoms: {len(viol)} ({len(vp)} primitive): {viol[:30]}", flush=True)
    json.dump({f"x_{i}": val[i] for i in range(NVARS)}, open('cand_chain.json', 'w'))
    print("wrote cand_chain.json", flush=True)

if __name__ == '__main__':
    main()
