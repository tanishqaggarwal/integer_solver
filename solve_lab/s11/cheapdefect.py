"""Minimise FAILING EQUATIONS, not bad checks: close expensive checks first and let
   1-equation checks absorb the damage."""
import sys, os, json, time, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
NAT = {u: len(L.var_atoms[u]) for u in range(L.NVARS)}
NEQ_OF = {a: len(L.atom2eq.get(a, {})) for a in range(L.NA)}


def load(name):
    v = [0] * L.NVARS
    for k, x in json.load(open(os.path.join(HERE, 'data', name))).items():
        v[int(k)] = int(x)
    fw.forward(v)
    return v


def cands(v, a, locked):
    out = []
    for u in L.avars[a]:
        if L.definer.get(u) is None and u not in locked and \
                not any(mm.count(u) > 1 for mm in L.polys[a]):
            out.append((u, None))
    out.sort(key=lambda kv: (NAT[kv[0]], kv[0]))
    try:
        hs, base = deep.handles(v, a, locked=locked)
        out += [(t, d) for t, d in sorted(hs, key=lambda kv: (NAT[kv[0]], kv[0]))]
    except Exception:
        pass
    return out


def try_close(v, a, locked):
    if fw.evalpoly(L.polys[a], v) == 0:
        return True
    for t, d in cands(v, a, locked):
        old = v[t]
        if d is None:
            x = fw.solve_lin(a, t, v)
            if x is None or x == old:
                continue
        else:
            bs = fw.evalpoly(L.polys[a], v)
            if not d or bs % d:
                continue
            x = old - bs // d
        v[t] = x
        fw.forward(v)
        if fw.evalpoly(L.polys[a], v) == 0:
            return True
        v[t] = old
        fw.forward(v)
    return False


def run(v, locked, rounds=25, verbose=True):
    best = (len(L.failing_eqs(L.all_atom_values(v))), [x for x in v])
    for rnd in range(rounds):
        bad = fw.bad_checks(v)
        if not bad:
            break
        # close the MOST EXPENSIVE first
        for a in sorted(bad, key=lambda a: -NEQ_OF[a]):
            try_close(v, a, locked)
        nb = fw.bad_checks(v)
        f = L.failing_eqs(L.all_atom_values(v))
        cost = sorted((NEQ_OF[a], a) for a in nb)
        if verbose:
            print(f"  rnd{rnd}: bad={len(nb)} failing={len(f)} score={L.NEQ-len(f)} costs={cost[:8]}", flush=True)
        if len(f) < best[0]:
            best = (len(f), [x for x in v])
        if not nb:
            break
    return best


if __name__ == '__main__':
    overall = None
    for name in ['closehit2.json', 'finish3.json', 'three.json', 'quad3_hit.json']:
        try:
            v = load(name)
        except Exception as e:
            print(name, 'skip', e)
            continue
        LOCK = {31339, 33708, 490, 91, 19750, 7497, 22820, 14853, 14393, 11436, 14515,
                16742, 22162, 30213, 8386, 21868, 16441, 28955, 2751, 18751}
        print(f"=== {name}: start failing={len(L.failing_eqs(L.all_atom_values(v)))}")
        b = run(v, LOCK)
        print(f"    -> best failing={b[0]} score={L.NEQ-b[0]}")
        if overall is None or b[0] < overall[0]:
            overall = b
    print(f"\nOVERALL best failing={overall[0]} score={L.NEQ-overall[0]}")
    json.dump({('x_%d' % i): overall[1][i] for i in range(L.NVARS)},
              open(os.path.join(HERE, 'data', 'cheapdefect_named.json'), 'w'))
