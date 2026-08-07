"""W stage 2: the EXACT integer congruence layer of every law block.
For each block: pair the c*N1 and c*N2 multiples, find the sum var, find the atom that
consumes it, and record liveness gate / outer multiplier / handle.  All from EQUATIONS.txt."""
import sys, os, re, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
import model
from collections import Counter
d = model.get(); A = d['atom_src']; AV = d['atom_vars']
byvar = {}
for i, vs in enumerate(AV):
    for v in vs: byvar.setdefault(v, []).append(i)
BL = json.load(open('w_blocks.json'))
blocks = BL['blocks']; PALIAS = set(BL['palias'])
def short(v): return [a for a in byvar.get(v, []) if len(A[a]) < 200]

shapes = Counter(); rows = []
nfail = 0
for b in blocks:
    m1 = {t[0]: t[1] for t in b['m1']}   # var -> coef  (coef * N1)
    m2 = {t[0]: t[1] for t in b['m2']}
    # sum atoms  x_Z - (x_Y1 + x_Y2)
    congs = []
    for y1 in m1:
        for a in short(y1):
            m = re.fullmatch(r'x_(\d+) - \(x_(\d+) \+ x_(\d+)\)', A[a])
            if not m: continue
            u, p1, p2 = int(m.group(1)), int(m.group(2)), int(m.group(3))
            other = p2 if p1 == y1 else (p1 if p2 == y1 else None)
            if other is None or other not in m2: continue
            congs.append({'Z': u, 'c1': m1[y1], 'c2': m2[other]})
    # a congruence may also be  c1*N1 alone summed with c2*N2 via the same pattern; also
    # the pair may be summed the other way round -- collected above symmetrically.
    b['congs'] = congs
    if len(congs) != 3: nfail += 1
    # now: what consumes each Z?
    for cg in congs:
        u = cg['Z']; cons = []
        for a in short(u):
            s = A[a]
            if re.fullmatch(r'x_\d+ - \(x_\d+ \+ x_\d+\)', s) and s.startswith('x_%d ' % u):
                continue
            cons.append((a, s))
        cg['cons'] = cons
        for a, s in cons:
            shapes[re.sub(r'x_\d+', 'V', re.sub(r'(?<![\dx_])\d{2,}', 'C', s))] += 1
print('blocks with != 3 congruences:', nfail)
print('congruence-consumer shapes:')
for k, v in shapes.most_common(20): print('  %5d  %s' % (v, k))
json.dump(blocks, open('w_blocks2.json', 'w'))
