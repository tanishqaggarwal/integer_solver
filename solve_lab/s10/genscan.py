"""S11 step 111: find EVERY zero-cost generator of the alpha-lattice.

lattice7 with nine measured generators finds that all 924 six-subsets and all 792
seven-subsets are solvable over Q and **none over Z**.  The obstruction is purely
integrality, and it comes from the two coarse directions:

    a0 moves only in multiples of 7376877   (x642, coupled to a6)
    a1 moves only in multiples of p         (x9413)

the fine generators for those being x7068 (cost 13) and x28730 (cost 16).  But those
were the only two ever tried.  a0 = x7068 - x2099 - 7376877*x642 also moves through
x2099, and a1 = x28730 - x9413*x17499 through x17499 -- so any free input feeding
those chains is a candidate generator, and the cost of each has to be MEASURED, since
the atom-count estimate was already wrong twice (x9118 and x8731 looked expensive and
are free).

Scan every free input in the structural support of the seven residual atoms: move it
by one, record the exact change to all seven alpha components, and count the failures
outside the twelve.

Usage: genscan.py START END
"""
import os, sys, time, json
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, ad
from frame2 import definer, ORDER, FREE, fwd
import suppfree
from chunk import sweep, load
P = ad.P
SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(base)
av0 = L.all_atom_values(base)
E = set(e for a in SEVEN for e in L.atom2eq[a])
BASE = L.NEQ - len(L.failing_eqs(av0))
print('witness %d; the twelve equations are %s' % (BASE, sorted(E)), flush=True)

_, freelist, SVS = suppfree.build(base, definer=definer, ORDER=ORDER, FREE=FREE,
                                  modp=None)
U = set()
for a in SEVEN:
    m = suppfree.atom_supp(a, base, SVS, modp=None)
    U |= {freelist[i] for i in range(len(freelist)) if (m >> i) & 1}
# a0 and a1 also move through x2099 and x17499, so add their supports
for t in (2099, 17499, 22665, 28599, 28961, 7075):
    m = SVS[t] if t < len(SVS) else 0
    U |= {freelist[i] for i in range(len(freelist)) if (m >> i) & 1}
U = sorted(U)
print('%d free inputs to scan' % len(U), flush=True)


def evaluate(u):
    v = list(base)
    v[u] = v[u] + 1
    fwd(v)
    aw = L.all_atom_values(v)
    d = [aw[a] - av0[a] for a in SEVEN]
    cost = len(set(L.failing_eqs(aw)) - E)
    return {'u': u, 'd': [str(x) for x in d], 'cost': cost,
            'nz': sum(1 for x in d if x)}


start = int(sys.argv[1]) if len(sys.argv) > 1 else 0
end = int(sys.argv[2]) if len(sys.argv) > 2 else len(U)
sweep('genscan', U, evaluate, start, min(end, len(U)), keyfn=str, budget=540)
rs = load('genscan')
free = [r for r in rs if r['cost'] == 0 and r['nz']]
print('\n%d scanned; %d move alpha at ZERO cost' % (len(rs), len(free)), flush=True)
NM = ['a22229', 'a22230', 'a35758', 'a35759', 'a35760', 'a35761', 'a35762']
for r in sorted(free, key=lambda r: r['u']):
    d = [int(x) for x in r['d']]
    print('   x%-6d  %s' % (r['u'], ', '.join(
        '%s%+d' % (NM[i], x) if abs(x) < 10 ** 12 else '%s%+.3g' % (NM[i], float(x))
        for i, x in enumerate(d) if x)), flush=True)
json.dump([{'u': r['u'], 'd': r['d']} for r in free],
          open(os.path.join(HERE, 'generators.json'), 'w'))
print('saved generators.json')
