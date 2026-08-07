"""W stage 15: SCREEN (not a proof).  The deliverable breaks the 2 off-pins of block E=7181
plus 5 handle atoms = 7 equations.  Which blocks touch the fewest equations through their
off-pins?  A lower count is a LEAD, not a price -- compensation inside an equation is exactly
what makes the real price smaller than the incidence."""
import sys, os, re, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
import model
from collections import Counter
d = model.get(); A = d['atom_src']; AV = d['atom_vars']; EQ = d['eq_terms']
byvar = {}
for i, vs in enumerate(AV):
    for v in vs: byvar.setdefault(v, []).append(i)
blocks = json.load(open('w_blocks4.json'))
eqs_of = {}
for i, (m, sq, tl) in enumerate(EQ):
    for c, a in tl: eqs_of.setdefault(a, set()).add(i)
def short(u): return [a for a in byvar.get(u, []) if len(A[a]) < 200]
PATS = [re.compile(r'x_(\d+) \* x_(\d+) - x_(\d+)$'), re.compile(r'x_(\d+) \* x_(\d+) \+ x_(\d+)$'),
        re.compile(r'(-?\d+) \* \(x_(\d+) \* x_(\d+)\) - x_(\d+)$'), re.compile(r'x_(\d+) \* x_(\d+) - (-?\d+) \* x_(\d+)$')]
rows = []
for b in blocks:
    L = b['L']
    NL = [int(re.match(r'x_(\d+)', A[a]).group(1)) for a in short(L)
          if re.fullmatch(r'x_\d+ - \(1 - x_%d\)' % L, A[a])][0]
    off = []
    for iv in (b['i5'], b['i6']):
        for a in short(iv):
            s = A[a]
            if not any(p.fullmatch(s) for p in PATS): continue
            if not ({NL, iv} <= set(int(x) for x in re.findall(r'x_(\d+)', s))): continue
            off.append(a)
    E = set()
    for a in off: E |= eqs_of.get(a, set())
    rows.append((len(E), b['E'], tuple(off)))
rows.sort()
print('equations touched by a block\'s two off-pins  (383 blocks):')
print('  histogram:', Counter(r[0] for r in rows).most_common(10))
print('  minimum 12:'); 
for r in rows[:12]: print('     %3d eqs   block E=%d   off-pins %s' % r)
me = [r for r in rows if r[1] == 7181]
print('  the deliverable\'s block E=7181:', me)
print('  rank of 7181 by incidence:', rows.index(me[0])+1, 'of', len(rows))
json.dump([{'neq': r[0], 'E': r[1], 'offpins': list(r[2])} for r in rows], open('w_price.json','w'))
