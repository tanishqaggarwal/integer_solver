"""Shared helpers for the quadrant re-solve.  Run from s9/ :  python3 quad/xxx.py"""
import pickle, sys, os, collections, json, time
S9 = '/home/user/integer_solver/solve_lab/s9'
os.chdir(S9)
sys.path.insert(0, S9)
import harness as H
exec(open('repair.py').read().split('if __name__')[0])

P = 2**256 - 2**32 - 977
NV = 38748
roots = pickle.load(open('roots.pkl', 'rb'))
atom2eq = pickle.load(open('atom2eq.pkl', 'rb'))
boolv = set(pickle.load(open('boolvars.pkl', 'rb')))
checks = [a for a in range(len(polys)) if a not in atom_out]
checkset = set(checks)
resid_poly = {a: (roots[a] if a in roots else polys[a]) for a in checks}
freeinp = [x for x in range(NV) if x not in definer]
freeset = set(freeinp)
bfree = [b for b in freeinp if b in boolv]

BEST = '../best/new_instance_partial_39022.json'


def nz_checks(v):
    """List of check atoms whose RESIDUAL (root for squares) is nonzero."""
    return sorted(a for a in checks if evalpoly(resid_poly[a], v) != 0)


def nz_gates(v):
    """Gate atoms that are nonzero (should be empty after a clean ripple)."""
    return sorted(a for a in atom_out if evalpoly(polys[a], v) != 0)


def eqs_of(atoms):
    s = set()
    for a in atoms:
        s.update(atom2eq.get(a, []))
    return sorted(s)


def fast_fails(v, nzatoms=None, codes=None):
    """Failing equations, evaluated only on equations touching a nonzero atom.
    Uses the RAW atom values (not roots) to decide which atoms are live."""
    if codes is None:
        codes = CODES
    live = [a for a in range(len(polys)) if evalpoly(polys[a], v) != 0]
    cand = eqs_of(live)
    return H.evaluate(codes, v, cand), live


_C = None
def CODESLOAD():
    global CODES
    CODES, _ = H.load_equations()
    return CODES
