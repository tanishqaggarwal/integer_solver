import time, json
from resources import marginal_window
out = {}
for mode in ('binary', 'wallace'):
    for w in (6, 7, 8, 9, 10):
        t0 = time.time()
        v, c, jr, mc = marginal_window(256, w, mode, neq=True, want_clique=True)
        out[f'{mode}_w{w}'] = dict(vars=v, couplers=c, jbits=jr, clique=mc)
        print(f"{mode:8s} w={w:2d}: vars={v:9,d} coupl={c:11,d} |J|=2^{jr} clique={mc} "
              f"({time.time()-t0:.0f}s)", flush=True)
        json.dump(out, open('window256_seq.json', 'w'), indent=1)
print("done")
