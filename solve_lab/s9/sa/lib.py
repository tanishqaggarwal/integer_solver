"""Shared library for the SA / stochastic bit search (agent 'sa').

Run everything with cwd = /home/user/integer_solver/solve_lab/s9 so that
repair.py / *.pkl resolve.
"""
import os, sys, pickle, json, time, random

S9 = '/home/user/integer_solver/solve_lab/s9'
if S9 not in sys.path:
    sys.path.insert(0, S9)
os.chdir(S9)

import harness as H

_ns = {}
exec(open(os.path.join(S9, 'repair.py')).read().split('if __name__')[0], _ns)
polys      = _ns['polys']
atom_out   = _ns['atom_out']
definer    = _ns['definer']
var_atoms  = _ns['var_atoms']
avars      = _ns['avars']
evalpoly   = _ns['evalpoly']
ripple     = _ns['ripple']
repair_loop= _ns['repair_loop']
solve_for  = _ns['solve_for']

NV = 38748
P  = 2**256 - 2**32 - 977
K1 = 33472904810391811973223207617762334363023286939839396241234196646906030803538671321618319
K2 = 42775533402728869434716629464193396056515231264222641773817154079369026410240838606908039

roots  = pickle.load(open('roots.pkl', 'rb'))
atom2eq= pickle.load(open('atom2eq.pkl', 'rb'))
boolv  = set(pickle.load(open('boolvars.pkl', 'rb')))
hits8599 = pickle.load(open('hits8599.pkl', 'rb'))

# residual polynomial per atom: use the degree-2 root for square check atoms
NA = len(polys)
resid_poly = [ (roots[a] if a in roots else polys[a]) for a in range(NA) ]

freeinp = [x for x in range(NV) if x not in definer]
freeset = set(freeinp)
bfree   = sorted(b for b in freeinp if b in boolv)   # the 1156 search-space bits

BASE_PATH = '../best/new_instance_partial_39022.json'
V0 = H.load_assignment(BASE_PATH)

CODES, VARSETS = H.load_equations()
NEQ = len(CODES)


# ---------------------------------------------------------------- fast objective
def nz_full(v):
    """All atoms whose residual (root for squares) is nonzero. O(all atoms)."""
    return set(a for a in range(NA) if evalpoly(resid_poly[a], v) != 0)


def nz_incremental(base_nz, changed_vars, v):
    """base_nz: nonzero set of the reference state; changed_vars: iterable of
    variables whose value differs from the reference state."""
    touched = set()
    for u in changed_vars:
        touched.update(var_atoms[u])
    nz = set(a for a in base_nz if a not in touched)
    for a in touched:
        if evalpoly(resid_poly[a], v) != 0:
            nz.add(a)
    return nz


def eqs_of(nz):
    s = set()
    for a in nz:
        s.update(atom2eq.get(a, ()))
    return s


def fast_fails(v, nz):
    """Exact evaluation of only the equations that can possibly fail."""
    ns = {'v': v, '__builtins__': {}}
    return [i for i in sorted(eqs_of(nz)) if eval(CODES[i], ns) != 0]


def true_fails(v):
    return H.evaluate(CODES, v)


# ---------------------------------------------------------------- decoder
def align(v):
    """The construct3 alignment tail (in place)."""
    ripple(v, {14853: v[12186]})
    ripple(v, {7068: v[2099] + 7376877 * v[642], 4432: v[19964] + v[28730]})
    ripple(v, {24548: v[25442]})
    return v


def decode(bits, alignment=True, k1=K1, base=None):
    """bits: iterable of boolean free inputs forced to 1.  Returns (v, changed_set)."""
    v = list(V0 if base is None else base)
    seeds = {b: 1 for b in bits}
    if alignment:
        seeds[5096] = k1
        seeds[33612] = 0
    ch = {}
    if seeds:
        c, _ = ripple(v, seeds)
        ch.update(c)
    if alignment:
        for s in ({14853: v[12186]},
                  {7068: v[2099] + 7376877 * v[642], 4432: v[19964] + v[28730]},
                  {24548: v[25442]}):
            c, _ = ripple(v, s)
            ch.update(c)
    return v, set(ch)


BASE_NZ = None   # filled by init_base()

def init_base():
    global BASE_NZ
    if BASE_NZ is None:
        BASE_NZ = nz_full(V0)
    return BASE_NZ


def score(bits, alignment=True, do_repair=0, k1=K1):
    """Return (nfail, failing_list, v, nz)."""
    init_base()
    v, ch = decode(bits, alignment=alignment, k1=k1)
    if do_repair:
        repair_loop(v, rounds=do_repair, verbose=False)
        nz = nz_full(v)
    else:
        nz = nz_incremental(BASE_NZ, ch, v)
    f = fast_fails(v, nz)
    return len(f), f, v, nz
