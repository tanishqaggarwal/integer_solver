"""Verify that model5's (scalar, linform-over-atoms) reproduces the raw equation value
exactly, on random integer assignments, for ALL equations."""
import pickle, random, sys, time
from bparse3 import parse_line, strip_g

W = '/home/user/integer_solver/solve_lab/agentB_work/'
M = pickle.load(open(W+'model5.pkl','rb'))
facs, atoms, eqs = M['facs'], M['atoms'], M['eqs']

def ev_ast(a, val):
    k = a[0]
    if k == 'g': return ev_ast(a[1], val)
    if k == 'n': return a[1]
    if k == 'v': return val[a[1]]
    if k == '+': return ev_ast(a[1][0], val) + ev_ast(a[1][1], val)
    if k == '-': return ev_ast(a[1][0], val) - ev_ast(a[1][1], val)
    if k == '*':
        r = 1
        for c in a[1]:
            r *= ev_ast(c, val)
            if r == 0: return 0
        return r
    raise ValueError(k)

def ev_poly(p, val):
    s = 0
    for m, c in p.items():
        t = c
        for v in m: t *= val[v]
        s += t
    return s

def main():
    rnd = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 12345)
    val = [rnd.randrange(-5, 6) for _ in range(38748)]
    faccache = {}
    bad = 0; t0 = time.time()
    for i, line in enumerate(open('/home/user/integer_solver/EQUATIONS.txt')):
        if not line.strip(): continue
        raw = ev_ast(parse_line(line), val)
        sc, L, k = eqs[i]
        if L is None:
            print("eq %d has no model" % i); bad += 1; continue
        s = 0
        for c, a in L:
            t = c
            for f in atoms[a]:
                fv = faccache.get(f)
                if fv is None:
                    fv = ev_poly(facs[f], val); faccache[f] = fv
                t *= fv
                if t == 0: break
            s += t
        pw = 1
        if k.startswith('pow'): pw = int(k[3:])
        elif k.startswith('same_pow'): pw = int(k[8:])
        mod = sc * (s ** pw)
        if mod != raw:
            bad += 1
            if bad < 6: print("MISMATCH eq %d: raw=%d model=%d kind=%s" % (i, raw, mod, k))
        if (i+1) % 10000 == 0:
            print("  %d checked, bad=%d %.1fs" % (i+1, bad, time.time()-t0), flush=True)
    print("TOTAL mismatches: %d / %d  (%.1fs)" % (bad, len(eqs), time.time()-t0))

if __name__ == '__main__':
    main()
