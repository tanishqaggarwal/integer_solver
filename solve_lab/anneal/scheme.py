#!/usr/bin/env python3
"""scheme.py -- given a machine, what is the best sound scheme and how many runs?

The deliverable the problem asks for is "n runs on hardware that exists".  This
turns the measured encoding costs into exactly that: for each machine, the most
bits of the scalar the annealer can own in one run, hence n.

Scheme (the only exact one -- see multirun/):  interval split.  The classical
outer loop fixes the top 256-mu bits of k, computes T' = T - k_hi*2^mu*G with one
EC subtraction, and hands the annealer a QUBO for the remaining mu bits.
    qubits per run = ceil(mu/w) * window_cost      runs n = 2^(256-mu)
"""
import json, math, os

# ---- hardware.  q = usable variables, clique = largest embeddable clique,
#      prec = usable bits of coupler precision.  Sources tracked in lit/SURVEY.md.
MACHINES = [
    # name                              quantum  q         clique  prec
    ("D-Wave Advantage (Pegasus P16)",  True,    5760,     180,    5),
    ("D-Wave Advantage2 (Zephyr Z15)",  True,    4400,     232,    5),
    ("D-Wave, hypothetical near-term",  True,    50000,    700,    6),
    ("Fujitsu Digital Annealer 3",      False,   100000,   100000, 64),
    ("Toshiba SQBM+ (classical)",       False,   1000000,  1000000, 32),
]


def load():
    enc = {}
    for f, tag in (('window256_seq.json', 'seq'), ('minvar.json', 'minvar')):
        if os.path.exists(f):
            for k, v in json.load(open(f)).items():
                enc[f"{tag}:{k}"] = v
    return enc


def capacity(machine, e):
    """largest number of LOGICAL variables of this encoding the machine can hold."""
    _, _, q, clique, prec = machine
    if e['jbits'] > prec:
        return 0, "coupler precision"
    chain = max(1.0, e['clique'] / 6.0) if clique < e['clique'] else 1.0
    if e['clique'] > clique:
        return 0, "clique too wide to embed"
    return int(q / chain), "ok"


if __name__ == '__main__':
    enc = load()
    if not enc:
        raise SystemExit("run remeasure.py / minvar.py first")
    print(f"{'machine':>34} {'best encoding':>22} {'logical cap':>12} "
          f"{'windows':>8} {'mu':>4} {'runs n':>10}")
    print("-" * 96)
    for m in MACHINES:
        best = None
        for name, e in enc.items():
            w = int(name.split('_w')[1].split('_')[0]) if '_w' in name else \
                int(name.split(':w')[1].split('_')[0])
            cap, why = capacity(m, e)
            nwin = cap // e['vars']
            mu = nwin * w
            if best is None or mu > best[0]:
                best = (mu, name, cap, nwin, why)
        mu, name, cap, nwin, why = best
        runs = f"2^{256-mu}" if mu < 256 else "1"
        note = "" if nwin else f"  <- {why}: cannot host one window"
        print(f"{m[0]:>34} {name.split(':')[-1]:>22} {cap:12,d} {nwin:8d} {mu:4d} {runs:>10}{note}")

    print()
    e = min(enc.values(), key=lambda v: v['vars'])
    w9 = enc.get('minvar:w9_uncompressed') or e
    print(f"smallest measured window: {min(v['vars'] for v in enc.values()):,d} logical variables")
    for target, label in ((128, "to beat Pollard rho (2^128)"), (256, "to solve in a single run")):
        need = math.ceil(target / 9) * min(v['vars'] for v in enc.values())
        print(f"  {label:<32}: mu >= {target:3d} bits -> >= {need:,d} logical variables")
