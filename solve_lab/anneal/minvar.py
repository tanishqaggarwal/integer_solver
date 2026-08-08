"""minvar.py -- the fewest-VARIABLES encoding of one comb window.

Chunk compression trades variables for smaller cliques and lower coupler
precision.  D-Wave wants small cliques; a fully-connected high-precision Ising
machine (Fujitsu Digital Annealer, Toshiba SQBM+) wants the opposite end of
that trade.  Measure both ends.
"""
import time, json
from resources import marginal_window

out = {}
for w in (7, 8, 9, 10):
    for ch, lbl in ((16, 'chunked'), (10**7, 'uncompressed')):
        t0 = time.time()
        v, c, jr, mc = marginal_window(256, w, 'binary', neq=True, want_clique=True, chunk=ch)
        out[f'w{w}_{lbl}'] = dict(vars=v, couplers=c, jbits=jr, clique=mc)
        print(f"w={w:2d} {lbl:>13}: vars={v:9,d} coupl={c:12,d} clique={mc:5d} |J|=2^{jr}"
              f"  ({time.time()-t0:.0f}s)", flush=True)
        json.dump(out, open('minvar.json', 'w'), indent=1)
print("done")
