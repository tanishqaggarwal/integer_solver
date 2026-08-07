"""Agent B: fast model-based scorer + residual analysis of an assignment."""
import pickle, json, sys, collections

W = '/home/user/integer_solver/solve_lab/agentB_work/'
_M = pickle.load(open(W+'model5.pkl','rb'))
facs, atoms, eqs = _M['facs'], _M['atoms'], _M['eqs']
NV = 38748

def load_assign(path):
    d = json.load(open(path))
    val = [0]*NV
    for k, v in d.items():
        val[int(k.split('_')[1])] = int(v)
    return val

def fac_values(val):
    fv = [0]*len(facs)
    for i, p in enumerate(facs):
        s = 0
        for m, c in p.items():
            t = c
            for v in m:
                t *= val[v]
                if t == 0: break
            s += t
        fv[i] = s
    return fv

def atom_values(fv):
    av = [0]*len(atoms)
    for i, a in enumerate(atoms):
        t = 1
        for f in a:
            t *= fv[f]
            if t == 0: break
        av[i] = t
    return av

def score(val):
    fv = fac_values(val); av = atom_values(fv)
    ok = 0; fail = []
    for e, (sc, L, k) in enumerate(eqs):
        s = 0
        for c, a in L: s += c*av[a]
        if s == 0: ok += 1
        else: fail.append(e)
    return ok, fail, fv, av

def main():
    path = sys.argv[1]
    val = load_assign(path)
    ok, fail, fv, av = score(val)
    print("model score: %d/%d  failing=%s" % (ok, len(eqs), fail[:20]))
    nzf = [i for i, x in enumerate(fv) if x]
    nza = [i for i, x in enumerate(av) if x]
    print("nonzero factors: %d   nonzero atoms: %d" % (len(nzf), len(nza)))
    # which nonzero atoms appear in failing equations
    inv = collections.Counter()
    for e in fail:
        for c, a in eqs[e][1]:
            if av[a]: inv[a] += 1
    print("nonzero atoms inside failing eqs:", len(inv))
    def s(p, lim=200):
        return ''.join('%+d*%s' % (c, '*'.join('x%d'%v for v in m) if m else '1')
                       for m, c in sorted(p.items()))[:lim]
    for a, n in inv.most_common(40):
        print("  atom %d (in %d failing eqs) value=%d" % (a, n, av[a]))
        for f in atoms[a]:
            print("      f%d = %s   -> %d" % (f, s(facs[f]), fv[f]))
    print("--- failing equations")
    for e in fail:
        L = eqs[e][1]
        tot = sum(c*av[a] for c, a in L)
        nz = [(c, a, av[a]) for c, a in L if av[a]]
        print("  eq %d |L|=%d residual=%d nonzero terms=%s" % (e, len(L), tot, nz))

if __name__ == '__main__':
    main()
