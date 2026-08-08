"""validate256.py -- build ONE comb window at the real field size to check the extrapolation."""
import time, json
from resources import marginal_window
out = {}
for w in (1, 8):
    t0 = time.time()
    v, c, jr = marginal_window(256, w, 'wallace')
    out[f'wallace_w{w}'] = dict(vars=v, couplers=c, jbits=jr, secs=round(time.time()-t0, 1))
    print(f"wallace s=256 w={w}: vars={v} couplers={c} |J|=2^{jr} ({time.time()-t0:.0f}s)", flush=True)
for w in (1, 8):
    t0 = time.time()
    v, c, jr = marginal_window(256, w, 'binary')
    out[f'binary_w{w}'] = dict(vars=v, couplers=c, jbits=jr, secs=round(time.time()-t0, 1))
    print(f"binary  s=256 w={w}: vars={v} couplers={c} |J|=2^{jr} ({time.time()-t0:.0f}s)", flush=True)
json.dump(out, open('window256.json', 'w'), indent=1)
