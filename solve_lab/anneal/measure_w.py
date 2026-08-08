import time, json, os
from resources import marginal_window
out = json.load(open('window256_neq.json')) if os.path.exists('window256_neq.json') else {}
for mode in ('binary', 'wallace'):
    for w in (4, 6, 7, 9, 11, 12):
        key = f'{mode}_w{w}'
        if key in out: continue
        t0 = time.time()
        v, c, jr = marginal_window(256, w, mode, neq=True)
        out[key] = dict(vars=v, couplers=c, jbits=jr)
        print(f"{mode:8s} w={w:2d}: vars={v:9,d} couplers={c:12,d} |J|=2^{jr} ({time.time()-t0:.0f}s)", flush=True)
        json.dump(out, open('window256_neq.json', 'w'), indent=1)
print("done")
