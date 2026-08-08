import time, json
from resources import marginal_window
out = {}
for mode in ('binary', 'wallace'):
    for w in (1, 8, 10):
        t0 = time.time()
        v, c, jr = marginal_window(256, w, mode, neq=True)
        out[f'{mode}_w{w}'] = dict(vars=v, couplers=c, jbits=jr)
        print(f"{mode:8s} s=256 w={w:2d} neq: vars={v:9,d} couplers={c:11,d} |J|=2^{jr} "
              f"({time.time()-t0:.0f}s)", flush=True)
json.dump(out, open('window256_neq.json', 'w'), indent=1)
