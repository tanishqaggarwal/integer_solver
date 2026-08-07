"""Channel measurement on the CORRECTED engine, at arbitrary configurations."""
import sys, os, json, re, collections, pickle
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine2 as E2, fast2

E_DIR = '/home/user/integer_solver/solve_lab/agentE_work'
F_DIR = '/home/user/integer_solver/solve_lab/agentF_work'
DELIV = '/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'
P = 115792089237316195423570985008687907853269984665640564039457584007908834671663
ROWS = [7389, 10187, 20212, 20215, 28647]     # E's cfg0 cluster
BAD8 = [23616, 23617, 36659, 36660, 36661, 36662, 36663, 36664]  # deliverable cluster


def load_vec(path=DELIV):
    d = json.load(open(path))
    v = [0] * H.NV
    for k, val in d.items():
        v[int(k.split('_')[1])] = int(val)
    return v


def isb(f):
    for i in H.occ[f]:
        t = re.sub(r'x_%d\b' % f, 'X', H.atoms[i])
        if t in ('X - X * X', 'X * X - X', 'X * (X - 1)', '2 * X * (1 - X)'):
            return True
    return False


_BOOLS = None
def bools():
    """The 256 boolean leaves (cone of ROWS), same knob set as before."""
    global _BOOLS
    if _BOOLS is None:
        import engine as E
        cand = sorted(set().union(*[set(E.cone(a)[1]) for a in ROWS]))
        _BOOLS = [f for f in cand if isb(f)]
    return _BOOLS


def measure(seed, atomset=None, targets=None):
    """Flip each boolean leaf 0->1 at this seed; return (v0, bad0, {leaf: signature})."""
    atomset = list(BAD8) if atomset is None else list(atomset)
    v0 = E2.forward(seed)
    bad0 = E2.badatoms(v0)
    B = targets if targets is not None else bools()
    sig = {}
    for f in B:
        if v0[f] != 0:
            sig[f] = 'ON'; continue
        b1, _ = fast2.resid_delta(v0, bad0, {f: 1})
        sig[f] = tuple((b1.get(a, 0) - bad0.get(a, 0)) % P for a in atomset)
    return v0, bad0, sig


def classes(sig):
    cls = collections.defaultdict(list)
    for f, c in sig.items():
        if c == 'ON':
            cls['ON'].append(f); continue
        if any(c):
            cls[c].append(f)
        else:
            cls['INERT'].append(f)
    return cls


def partition(sig):
    """Just the block structure (list of frozensets), ON/INERT included as blocks."""
    cls = classes(sig)
    return sorted((frozenset(v) for v in cls.values()), key=lambda s: (-len(s), min(s)))


def tree():
    t = json.load(open(os.path.join(F_DIR, 'tree96.json')))
    return t, {k: set(v['gsup']) for k, v in t.items()}
