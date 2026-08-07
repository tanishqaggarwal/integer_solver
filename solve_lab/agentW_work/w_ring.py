"""W stage 3: THE RING.  For every one of the 1149 congruence atoms, decompose exactly:
      sgn * ( a * (L * Z) )  =  c * H ,   Z = c1*N1 + c2*N2
and resolve H.  Determines whether the gadget condition is  P | .. ,  c*P | .. , or weaker."""
import sys, os, re, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
import model
from collections import Counter
d = model.get(); A = d['atom_src']; AV = d['atom_vars']
byvar = {}
for i, vs in enumerate(AV):
    for v in vs: byvar.setdefault(v, []).append(i)
blocks = json.load(open('w_blocks2.json'))
PALIAS = set(json.load(open('w_blocks.json'))['palias'])
def short(v): return [a for a in byvar.get(v, []) if len(A[a]) < 200]

PATS = [
 (re.compile(r'x_(\d+) \* x_(\d+) - x_(\d+)$'),        lambda m: (1, int(m.group(1)), int(m.group(2)), 1, int(m.group(3)))),
 (re.compile(r'x_(\d+) \* x_(\d+) \+ x_(\d+)$'),       lambda m: (1, int(m.group(1)), int(m.group(2)), -1, int(m.group(3)))),
 (re.compile(r'(-?\d+) \* \(x_(\d+) \* x_(\d+)\) - x_(\d+)$'), lambda m: (int(m.group(1)), int(m.group(2)), int(m.group(3)), 1, int(m.group(4)))),
 (re.compile(r'x_(\d+) \* x_(\d+) - (-?\d+) \* x_(\d+)$'), lambda m: (1, int(m.group(1)), int(m.group(2)), int(m.group(3)), int(m.group(4)))),
]
gatecnt = Counter(); hcls = Counter(); mult_out = Counter(); mult_h = Counter()
recs = []; unresolved = []
for b in blocks:
    for cg in b['congs']:
        Z = cg['Z']; got = None
        for a, s in cg['cons']:
            for pat, f in PATS:
                m = pat.fullmatch(s)
                if m:
                    aout, v1, v2, chand, H = f(m)
                    L = v1 if v2 == Z else (v2 if v1 == Z else None)
                    if L is None: continue
                    got = dict(atom=a, aout=aout, L=L, Z=Z, chand=chand, H=H,
                               c1=cg['c1'], c2=cg['c2'])
            if got: break
        if not got: unresolved.append((b['E'], Z)); continue
        # resolve H:  expect  x_H - x_Palias * x_u   (u free handle)
        hd = None
        for a in short(got['H']):
            m = re.fullmatch(r'x_%d - x_(\d+) \* x_(\d+)' % got['H'], A[a])
            if m:
                p, q = int(m.group(1)), int(m.group(2))
                if p in PALIAS: hd = ('P*u', q)
                elif q in PALIAS: hd = ('P*u', p)
                else: hd = ('prod', (p, q))
            elif re.fullmatch(r'x_%d - .*' % got['H'], A[a]): hd = hd or ('other', A[a])
        got['H_kind'] = hd[0] if hd else 'UNDEF'
        got['u'] = hd[1] if hd and hd[0] == 'P*u' else None
        hcls[got['H_kind']] += 1
        mult_out[got['aout']] += 1; mult_h[got['chand']] += 1
        gatecnt[got['L']] += 0
        recs.append(got)
        cg['ring'] = got
print('congruence atoms resolved:', len(recs), ' unresolved:', len(unresolved))
print('handle kind:', dict(hcls))
print('outer multipliers a (on L*Z):', mult_out.most_common(6))
print('handle multipliers c:', mult_h.most_common(6))
# is the handle u private (appears only in its own defining atom + big eq atoms)?
priv = Counter()
for r in recs:
    if r['u'] is None: continue
    defs = [a for a in short(r['u']) if re.fullmatch(r'x_\d+ - .*', A[a])]
    priv[len(byvar[r['u']]), len(defs)] += 1
print('handle u  (total atoms, short-defining atoms) histogram:', priv.most_common(8))
json.dump(blocks, open('w_blocks3.json', 'w'))
