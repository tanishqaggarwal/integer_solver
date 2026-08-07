"""W stage 8: the OFF-PINS.  Each block also carries  a*( (1-L) * i5 ) = c*P*u  and the same
for i6.  So gate=0 does not free the output -- it pins it to 0 mod P."""
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
PALIAS = set(json.load(open('w_blocks.json'))['palias'])
def short(v): return [a for a in byvar.get(v, []) if len(A[a]) < 200]
PATS = [
 (re.compile(r'x_(\d+) \* x_(\d+) - x_(\d+)$'),                lambda m: (1, m.group(1), m.group(2), 1, m.group(3))),
 (re.compile(r'x_(\d+) \* x_(\d+) \+ x_(\d+)$'),               lambda m: (1, m.group(1), m.group(2), -1, m.group(3))),
 (re.compile(r'(-?\d+) \* \(x_(\d+) \* x_(\d+)\) - x_(\d+)$'), lambda m: (m.group(1), m.group(2), m.group(3), 1, m.group(4))),
 (re.compile(r'x_(\d+) \* x_(\d+) - (-?\d+) \* x_(\d+)$'),     lambda m: (1, m.group(1), m.group(2), m.group(3), m.group(4))),
]
res = Counter(); bad = []
for b in blocks:
    L = b['L']
    # NOT-gate var:  x_nl - (1 - x_L)
    nl = [int(re.match(r'x_(\d+)', A[a]).group(1)) for a in short(L)
          if re.fullmatch(r'x_\d+ - \(1 - x_%d\)' % L, A[a])]
    if len(nl) != 1: res['no-NOT'] += 1; bad.append((b['E'], 'notgate', nl)); continue
    NL = nl[0]
    for slot, iv in (('i5', b['i5']), ('i6', b['i6'])):
        found = []
        for a in short(iv):
            s = A[a]
            for pat, f in PATS:
                m = pat.fullmatch(s)
                if not m: continue
                aout, v1, v2, ch, H = f(m); v1, v2, H = int(v1), int(v2), int(H)
                if {v1, v2} != {NL, iv}: continue
                # H must be P * free
                hk = None
                for a2 in short(H):
                    m2 = re.fullmatch(r'x_%d - x_(\d+) \* x_(\d+)' % H, A[a2])
                    if m2:
                        p, q = int(m2.group(1)), int(m2.group(2))
                        if p in PALIAS: hk = ('P*u', q)
                        elif q in PALIAS: hk = ('P*u', p)
                        else: hk = ('prod', None)
                found.append((int(aout), int(ch), hk[0] if hk else 'UNDEF', hk[1] if hk else None))
        if len(found) == 1 and found[0][2] == 'P*u':
            u = found[0][3]
            res['ok, handle private' if len(byvar[u]) <= 6 else 'ok'] += 1
        else:
            res['BAD'] += 1; bad.append((b['E'], slot, iv, found))
print('off-pin census over 766 output slots:', res.most_common())
print('bad:', len(bad), bad[:5])
