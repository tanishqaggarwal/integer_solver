#!/usr/bin/env python3
"""measure_grid.py -- MEASURED per-window cost across field sizes, for the runs table."""
import sys, json, time, os
sys.path.insert(0, '.')
from resources import marginal_window

OUT = 'synth/window_grid.json'
grid = json.load(open(OUT)) if os.path.exists(OUT) else {}
for s in (16, 24, 32, 48, 64, 96, 128):
    for w in (1, 2, 4, 6, 8):
        key = f'{s}_{w}'
        if key in grid: continue
        t = time.time()
        try:
            v, c, jr, mc = marginal_window(s, w, 'binary', neq=True, want_clique=True, chunk=16)
            grid[key] = dict(s=s, w=w, vars=v, couplers=c, jbits=jr, clique=mc)
            print(f"s={s:3d} w={w}: {v:7d} vars clique={mc} |J|=2^{jr} ({time.time()-t:.0f}s)", flush=True)
            json.dump(grid, open(OUT, 'w'), indent=1)
        except Exception as e:
            print(f"s={s} w={w}: {e}")
print("done")
