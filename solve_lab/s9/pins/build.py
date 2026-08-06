"""Canonical builder over the 256-bit message + pin closure, and defect measurement."""
import pickle, collections, sys, os, json, time
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s9')
os.chdir('/home/user/integer_solver/solve_lab/s9')
import harness as H
exec(open('repair.py').read().split('if __name__')[0])

P = 2**256 - 2**32 - 977
NV = 38748
roots = pickle.load(open('roots.pkl', 'rb'))
checks = [a for a in range(len(polys)) if a not in atom_out]
rp = {a: (roots[a] if a in roots else polys[a]) for a in checks}
pins = pickle.load(open('pins/pins.pkl', 'rb'))
bitpins = collections.defaultdict(list)
for p_ in pins: bitpins[p_['G']].append(p_)
BITS = sorted(bitpins)
freeset = set(x for x in range(NV) if x not in definer)
K1 = 33472904810391811973223207617762334363023286939839396241234196646906030803538671321618319
K2 = 42775533402728869434716629464193396056515231264222641773817154079369026410240838606908039

BASE = H.load_assignment('../best/new_instance_partial_39022.json')

def nz(v):
    return sorted(a for a, Pp in rp.items() if evalpoly(Pp, v) != 0)

def close_mirrors(v):
    ripple(v, {7068: v[2099] + 7376877*v[642], 4432: v[19964] + v[28730]})
    ripple(v, {24548: v[25442], 14853: v[1308]})

def build(bitset, extra=None, close=True):
    """bitset: iterable of bit vars to set to 1 (all others 0)."""
    v = list(BASE)
    s = {}
    on = set(bitset)
    for b in BITS:
        s[b] = 1 if b in on else 0
    for b in BITS:
        for pn in bitpins[b]:
            s[pn['h']] = 0
            if b in on: s[pn['B']] = pn['HUGE']
    if extra: s.update(extra)
    ripple(v, s)
    if close: close_mirrors(v)
    return v

def defects(v):
    A = evalpoly(polys[22229], v)
    B = evalpoly(polys[22231], v)
    return A % P, B % P

if __name__ == '__main__':
    t0 = time.time()
    v = build([2081, 24601])
    n = nz(v)
    print('base bitset {2081,24601}: residual atoms', n)
    print('defects mod p:', defects(v))
    codes, _ = H.load_equations()
    f = H.evaluate(codes, v)
    print(f'EQUATIONS satisfied {len(codes)-len(f)}/{len(codes)} ({len(f)} failing) {f}')
    print(f'{time.time()-t0:.0f}s')
