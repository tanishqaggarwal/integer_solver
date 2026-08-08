#!/usr/bin/env python3
"""measure.py -- marginal cost of one window at s = 256, one config per process.

Run as   python3 measure.py run <json-key> <kv,kv,...>   for a single config
(so the interpreter's memory is released between configs), or
         python3 measure.py all
to drive the whole grid through subprocesses and collect arith/win256.json.
"""
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'win256.json')


def parse(kv):
    d = {}
    for it in kv.split(','):
        if not it:
            continue
        k, v = it.split('=')
        d[k] = (v == 'True') if v in ('True', 'False') else (
            v if k in ('mode', 'onehot', 'kind') else int(v))
    return d


def run_one(cfg):
    from marginal import window, semaev_step
    kind = cfg.pop('kind', 'window')
    fn = semaev_step if kind == 'semaev' else window
    st = fn(256, **cfg)
    return dict(st)


GRID = []
for w in (6, 8, 9, 10, 11):
    GRID.append((f"base_w{w}", dict(w=w, mux=False, onehot='square', kdepth=0)))
for w in (6, 8, 9, 10):
    GRID.append((f"mux_w{w}", dict(w=w, mux=True, onehot='square', kdepth=0)))
for w in (8, 10, 11, 12):
    GRID.append((f"muxtree_w{w}", dict(w=w, mux=True, onehot='tree', kdepth=0)))
for w in (8, 10, 11, 12):
    GRID.append((f"muxtreekara_w{w}",
                 dict(w=w, mux=True, onehot='tree', kdepth=4)))
for w in (9, 11, 12, 13):
    GRID.append((f"full_w{w}", dict(w=w, mux=True, onehot='tree', kdepth=4,
                                    signed=True)))
for ch in (32, 64):
    GRID.append((f"full_w12_c{ch}", dict(w=12, mux=True, onehot='tree',
                                         kdepth=4, signed=True, chunk=ch)))
for w in (11, 13):
    GRID.append((f"semaev_w{w}", dict(kind='semaev', w=w, mux=True,
                                      onehot='tree', kdepth=4)))
GRID.append(("basekara_w9", dict(w=9, mux=False, onehot='square', kdepth=4)))


if __name__ == '__main__':
    if sys.argv[1] == 'run':
        sys.path.insert(0, HERE)
        print(json.dumps(run_one(parse(sys.argv[2]))))
    else:
        out = json.load(open(OUT)) if os.path.exists(OUT) else {}
        for key, cfg in GRID:
            if key in out:
                continue
            kv = ','.join(f"{k}={v}" for k, v in cfg.items())
            t0 = time.time()
            r = subprocess.run([sys.executable, __file__, 'run', kv],
                               capture_output=True, text=True, cwd=HERE)
            if r.returncode != 0:
                print(f"{key:22s} FAILED: {r.stderr.strip().splitlines()[-1:]}",
                      flush=True)
                continue
            st = json.loads(r.stdout.strip().splitlines()[-1])
            st['cfg'] = kv
            out[key] = st
            print(f"{key:22s} vars={st['vars']:10,d} coupl={st['couplers']:12,d} "
                  f"AND={st['and_vars']:9,d} |J|=2^{st['dynamic_range_bits']:<3d}"
                  f" ({time.time()-t0:.0f}s)", flush=True)
            json.dump(out, open(OUT, 'w'), indent=1)
        print("done")
