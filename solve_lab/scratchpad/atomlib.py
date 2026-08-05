"""Shared atom/equation library for the integer_solver lab.
Load poly atoms (each is list of [varlist,coeff]), evaluate at assignments,
find nonzero atoms, trace variable->atom incidence."""
import json, re
p = 2**256 - 2**32 - 977
NVARS = 38748
BASE = '/home/user/integer_solver/solve_lab'

# --- load poly atoms ---
ATOMS = []          # list of poly (list of [varlist, coeff])
ATOM_REPR = []
ATOM_EQS = []       # list of eqs each atom appears in
with open(BASE + '/atoms/poly_atoms.jsonl') as f:
    for line in f:
        d = json.loads(line)
        ATOMS.append(d['poly'])
        ATOM_REPR.append(d.get('repr', ''))
        ATOM_EQS.append(d.get('eqs', []))
NATOM = len(ATOMS)

# atom -> set of vars
ATOM_VARS = []
for poly in ATOMS:
    s = set()
    for varlist, coeff in poly:
        for v in varlist:
            s.add(v)
    ATOM_VARS.append(s)

# var -> list of atom indices
from collections import defaultdict
VAR_ATOMS = defaultdict(list)
for ai, s in enumerate(ATOM_VARS):
    for v in s:
        VAR_ATOMS[v].append(ai)

def eval_atom(ai, v):
    tot = 0
    for varlist, coeff in ATOMS[ai]:
        term = coeff
        for var in varlist:
            term *= v[var]
        tot += term
    return tot

def nonzero_atoms(v, mod=None):
    out = []
    for ai in range(NATOM):
        val = eval_atom(ai, v)
        if mod is not None:
            val %= mod
        if val != 0:
            out.append((ai, val))
    return out

def load_json(path):
    if not path.startswith('/'):
        path = BASE + '/' + path
    d = json.load(open(path))
    v = [0] * NVARS
    for k, val in d.items():
        idx = int(k[2:]) if k.startswith('x_') else int(k)
        v[idx] = int(val)
    return v

# --- equation loader ---
VAR_RE = re.compile(r'x_(\d+)')
_eqcode = None
_eqvars = None
def load_eqs():
    global _eqcode, _eqvars
    if _eqcode is not None:
        return _eqcode, _eqvars
    lines = [L for L in open(BASE + '/../EQUATIONS.txt').read().split('\n') if L.strip()]
    _eqcode = [compile(VAR_RE.sub(r'v[\1]', L.rsplit('=', 1)[0]), '<e>', 'eval') for L in lines]
    _eqvars = [frozenset(int(m) for m in VAR_RE.findall(L)) for L in lines]
    return _eqcode, _eqvars

def eval_eqs(v, idxs=None):
    code, _ = load_eqs()
    ns = {'v': v, '__builtins__': {}}
    if idxs is None:
        idxs = range(len(code))
    return {i: eval(code[i], ns) for i in idxs}

def failing(v):
    code, _ = load_eqs()
    ns = {'v': v, '__builtins__': {}}
    return [i for i in range(len(code)) if eval(code[i], ns) != 0]

if __name__ == '__main__':
    import sys
    v = load_json(sys.argv[1])
    nz = nonzero_atoms(v)
    print(f"{len(nz)} nonzero atoms at {sys.argv[1]}")
    for ai, val in nz:
        vp = val % p
        print(f"  atom {ai}: {ATOM_REPR[ai]!r}  val%p={vp}  ({'==0 mod p' if vp==0 else 'NONZERO mod p'})  vars={sorted(ATOM_VARS[ai])}  in {len(ATOM_EQS[ai])} eqs")
