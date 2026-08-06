"""All-zero branch: the forced OR gate can ONLY be satisfied by one of the 256 gate-bits,
and every gate-bit carries exactly two load pins.  So the branch reduces to a 256-way choice:
which gate-bit's pins are cheapest to close?  Scan all of them, both with pins left open and
with pins closed, and keep the best."""
import pickle, sys, time
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
P = 2**256 - 2**32 - 977
d = pickle.load(open('atoms.pkl', 'rb')); src = d['atom_src']
a2e = {a: set(e) for a, e in pickle.load(open('atom2eq.pkl', 'rb')).items()}
roots = pickle.load(open('roots.pkl', 'rb'))
pins = pickle.load(open('pins/pins.pkl', 'rb'))
bypin = {}
for pn in pins:
    bypin.setdefault(pn['G'], []).append(pn)
GATEBITS = sorted(bypin)
CODES, _ = H.load_equations()


def allnz(v):
    return sorted(a for a, Pp in enumerate(polys) if evalpoly(Pp, v) != 0)


def fast_fail(v):
    """equations that can fail = those containing a nonzero atom; evaluate only those, exactly"""
    nz = allnz(v)
    eqs = set()
    for a in nz: eqs |= a2e.get(a, set())
    ns = {'v': v, '__builtins__': {}}
    return [i for i in sorted(eqs) if eval(CODES[i], ns) != 0]


if __name__ == '__main__':
    v0 = H.load_assignment('zero_out.json')
    ripple(v0, {22162: 0, 30213: 0})
    # validate the fast objective once against the real evaluator
    a = len(fast_fail(v0)); b = len(H.evaluate(CODES, v0))
    print(f'fast objective validation: {a} vs real {b}  -> {"OK" if a == b else "MISMATCH"}')
    assert a == b
    best = None
    t0 = time.time()
    for i, G in enumerate(GATEBITS):
        for close in (False, True):
            w = list(v0)
            seeds = {G: 1}
            if close:
                for pn in bypin[G]: seeds[pn['B']] = pn['HUGE']
            ripple(w, seeds)
            repair_loop(w, rounds=10, verbose=False)
            f = fast_fail(w)
            if best is None or len(f) < best[0]:
                best = (len(f), G, close, allnz(w))
                print(f'   new best: x_{G} close_pins={close} -> {len(f)} failing, atoms {best[3]}', flush=True)
        if i % 50 == 0: print(f'   ...{i}/{len(GATEBITS)}  {time.time()-t0:.0f}s', flush=True)
    n, G, close, nz = best
    print(f'\nBEST gate-bit x_{G} (close_pins={close}): {n} failing, atoms {nz}')
    w = list(v0); seeds = {G: 1}
    if close:
        for pn in bypin[G]: seeds[pn['B']] = pn['HUGE']
    ripple(w, seeds); repair_loop(w, rounds=10, verbose=False)
    print('verify with real checker:', len(H.evaluate(CODES, w)), 'failing')
    for a in allnz(w):
        print(f'    atom {a} ({len(a2e.get(a,[]))} eqs) gate={atom_out.get(a)}: {src[a][:95]}')
    H.save_assignment(w, 'zero5_out.json')
