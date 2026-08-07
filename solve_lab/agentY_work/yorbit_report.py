#!/usr/bin/env python3
"""Build the per-target orbit table for RESUME_Y.md section 6.1 and splice it in.
A target counts as DONE only with a DONE line for every one of sizes 2, 3, 4, and the
candidate counts must equal C(256,b) exactly."""
import os, re, glob
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__))
NAMES = ['negT','lamT','neglamT','lam2T','neglam2T',
         'c_negT','c_lamT','c_neglamT','c_lam2T','c_neglam2T']
SCAL = {'negT':'-k','lamT':'lam*k','neglamT':'-lam*k','lam2T':'lam^2*k','neglam2T':'-lam^2*k',
        'c_negT':'(2^256-1)+k','c_lamT':'(2^256-1)-lam*k','c_neglamT':'(2^256-1)+lam*k',
        'c_lam2T':'(2^256-1)-lam^2*k','c_neglam2T':'(2^256-1)+lam^2*k'}
want = {b: comb(256, b) for b in (2, 3, 4)}
rows, ndone = [], 0
for nm in NAMES:
    f = os.path.join(HERE, 'rep_orbit_%s.txt' % nm)
    got, hits, zero, secs = {}, 0, 0, 0.0
    if os.path.exists(f):
        for l in open(f):
            m = re.match(r'DONE size=(\d+) range=\[0,256\) n=(\d+) zero=(\d+) ([\d.]+)s', l)
            if m:
                got[int(m.group(1))] = int(m.group(2)); zero += int(m.group(3)); secs += float(m.group(4))
            if l.startswith('HIT'): hits += 1
    ok = all(b in got and got[b] == want[b] for b in (2, 3, 4))
    ndone += ok
    if ok:
        rows.append('| `%s` | `%s` | **exhausted** | %d | %d | %.1f s |' % (nm, SCAL[nm], hits, zero, secs))
    else:
        have = ','.join(str(b) for b in sorted(got)) or 'none'
        rows.append('| `%s` | `%s` | *incomplete* — sizes done: %s | %d | %d | — |' % (nm, SCAL[nm], have, hits, zero))
tbl = ['| target | scalar | `|S| <= 8` | hits | degenerate | time |',
       '|---|---|---|---|---|---|'] + rows
tbl.append('')
tbl.append('**%d of 10 targets complete.** Each complete target scanned %d candidates '
           '(`C(256,2)+C(256,3)+C(256,4)`), counts checked exactly against the binomials; '
           'combined with the `|S| <= 4` probe of §5 this exhausts `|S| <= 8` on that target. '
           'Targets marked *incomplete* are **not** exhausted and must not be quoted as such.'
           % (ndone, sum(want.values())))
block = '\n'.join(tbl)
p = os.path.join(HERE, 'RESUME_Y.md')
src = open(p).read()
new = re.sub(r'<!--ORBIT-->.*?(?=\n### 6\.2)', '<!--ORBIT-->\n' + block + '\n', src, flags=re.S)
open(p, 'w').write(new)
print(block)
print('\nspliced into RESUME_Y.md (%d/10 complete)' % ndone)
