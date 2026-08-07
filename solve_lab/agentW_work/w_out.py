"""W stage 7: ALL consumers of the output pair (i5,i6) of every block.  If the only places
they occur are the block's own law atoms and the two gated mux products L*i5, L*i6, then
gate=0 kills both the law AND the output -> 'gate off' is not an exploitable family."""
import sys, os, re, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
import model
from collections import Counter
d = model.get(); A = d['atom_src']; AV = d['atom_vars']
byvar = {}
for i, vs in enumerate(AV):
    for v in vs: byvar.setdefault(v, []).append(i)
blocks = json.load(open('w_blocks4.json'))
bad = []; census = Counter(); nlaw = Counter()
for b in blocks:
    L = b['L']
    for slot, iv in (('i5', b['i5']), ('i6', b['i6'])):
        others = []
        for a in byvar[iv]:
            s = A[a]
            if len(s) >= 250: continue                       # equation-level composites
            if re.fullmatch(r'x_\d+ - x_%d \* x_%d' % (L, iv), s) or \
               re.fullmatch(r'x_\d+ - x_%d \* x_%d' % (iv, L), s):
                census['gated mux product'] += 1; continue
            # the block's own law atoms
            if slot == 'i5' and re.fullmatch(r'x_%d - \(x_\d+ \+ x_%d\)|x_%d - \(x_%d \+ x_\d+\)'
                                             % (b['D'], iv, b['D'], iv), s):
                census['law: D sum'] += 1; continue
            if slot == 'i5' and re.fullmatch(r'x_%d - \(x_\d+ - x_%d\)' % (b['Jx'] if False else 0, iv), s):
                pass
            others.append((a, s))
        # classify the leftovers
        for a, s in others:
            tag = 'LAW' if any(('x_%d' % w) in s for w in (b['A'], b['B'], b['G'], b['E'], b['D'], b['N1'], b['N2'])) \
                  or re.fullmatch(r'x_\d+ - \(x_\d+ [-+] x_%d\)' % iv, s) \
                  or re.fullmatch(r'x_\d+ - \(x_%d [-+] x_\d+\)' % iv, s) else 'OTHER'
            census[tag] += 1
            if tag == 'OTHER': bad.append((b['E'], slot, iv, a, s[:140]))
print('consumer census of the 383 output pairs:', census.most_common())
print('OTHER consumers (would break the argument):', len(bad))
for t in bad[:12]: print('  ', t)
# and: does anything OUTSIDE the block read i5/i6 raw?  list the distinct shapes of every
# short atom touching i5/i6 that is not the gated product
shapes = Counter()
for b in blocks:
    for iv in (b['i5'], b['i6']):
        for a in byvar[iv]:
            s = A[a]
            if len(s) >= 250: continue
            shapes[re.sub(r'x_\d+', 'V', re.sub(r'(?<![\dx_])\d{2,}', 'C', s))] += 1
print('short-atom shapes touching i5/i6:', shapes.most_common(10))
