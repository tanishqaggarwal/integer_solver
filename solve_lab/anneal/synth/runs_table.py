#!/usr/bin/env python3
"""runs_table.py -- runs = 2^(bits - mu), mu maximized per machine, from MEASURED costs.

per-window cost measured directly at each field size (synth/window_grid.json for
s<=128, window256_seq.json for s=256). A comb resolving mu scalar bits needs
ceil(mu/w) windows, each doing full s-bit field arithmetic (s = key bit-size).
Reducing runs === raising mu === shrinking one window.  Baseline encoder here;
the synth/mincost track improves per-window ~3x, which shifts every mu up.
"""
import json, math

G = json.load(open('synth/window_grid.json'))
W256 = json.load(open('window256_seq.json'))

def windows(s):
    if s == 256:
        return {int(k.split('_w')[1]): (v['vars'], v['clique'], v['jbits'])
                for k, v in W256.items() if k.startswith('binary_w')}
    return {v['w']: (v['vars'], v['clique'], v['jbits'])
            for k, v in G.items() if v['s'] == s}

def mu_max(budget, s, prec):
    ws = windows(s)
    best = 0, None
    for w, (pw, cl, jb) in ws.items():
        if jb > prec: continue                      # coupler-precision gate
        M = int(budget // pw)
        if M < 1: continue
        mu = min(w * M, s)
        if mu > best[0]: best = mu, w
    return best

MACHINES = [
    ("D-Wave Advantage2 (Zephyr)",   4400,   5),
    ("D-Wave Advantage  (Pegasus)",  5760,   5),
    ("near-term D-Wave (~50k)",      50000,  6),
    ("Fujitsu DA3 (8k, full, 64b)",  8192,   64),
    ("Toshiba SB (1e5, full, fp)",   100000, 32),
    ("idealized 1e6-qubit",          1000000, 64),
]

if __name__ == '__main__':
    print("RUNS to recover a b-bit ECDLP key = 2^(b - mu), mu = scalar bits annealed per run")
    print("(measured baseline encoder; mu capped by qubit budget AND coupler precision)\n")
    for bits in (32, 48, 64, 128, 256):
        print(f"=== {bits}-bit key (field s={bits}, one window = "
              f"{min(v[0] for v in windows(bits).values()):,d} qubits min) ===")
        print(f"  {'machine':>32} {'qubits':>8} {'prec':>5} {'mu/run':>7} {'runs':>12}")
        for name, q, prec in MACHINES:
            mu, w = mu_max(q, bits, prec)
            if mu >= bits:  tag = "1 (whole key)"
            elif mu == 0:   tag = f"2^{bits} (no fit)"
            else:           tag = f"2^{bits-mu}"
            print(f"  {name:>32} {q:8d} {prec:5d} {mu:7d} {tag:>12}")
        print()
