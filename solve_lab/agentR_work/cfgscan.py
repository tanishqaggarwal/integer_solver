#!/usr/bin/env python3
"""Score real full assignments for chosen selector configurations.

For a configuration (a set of conditional-pin booleans turned on, plus its tree partner flags),
build a complete assignment with agent F's chain-repair `gs2.solve` (imported read-only), then
score it with agent F's exact evaluator.  Checkpointed to runs/cfgscan.json.

This asks the ONE question that can beat the deliverable without solving the target equation:
does the number of failing equations depend on the configuration, and can it drop below 7?
"""
import sys, os, json, time, pickle, itertools, random
FW = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, FW)
import gs2
from fwd import NV
E = gs2.E
T1 = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
T2 = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
pins = json.load(open(os.path.join(FW, 'pins.json')))
sup = pickle.load(open(os.path.join(FW, 'supp.pkl'), 'rb'))
TREE = {}
for b in pins:
    if int(b) in sup['7715']: TREE[int(b)] = 'A'
    elif int(b) in sup['34554']: TREE[int(b)] = 'B'
PARTNER = {'A': 5090, 'B': 22106}
PIN = {22162: T1, 30213: T2, 24468: T1, 18956: T2}

def run_cfg(bits):
    v = [0] * NV
    fr = set(PIN) | set(bits)
    for k, x in PIN.items(): v[k] = x
    for b in bits: v[b] = 1
    for b in bits:
        pa = PARTNER[TREE[b]]; v[pa] = 1; fr.add(pa)
    v, ok = gs2.solve(v, verbose=False, frozen=set(fr))
    r = E.run(v); bad = E.score(r)
    return 39033 - len(bad), sum(1 for x in r if x), ok, v

OUT = 'runs/cfgscan.json'
res = json.load(open(OUT)) if os.path.exists(OUT) else {}

def key(bits): return ','.join(map(str, sorted(bits)))

def do(bits, tag=''):
    k = key(bits)
    if k in res: return res[k]
    t = time.time()
    try:
        sc, nz, ok, v = run_cfg(bits)
    except Exception as e:
        res[k] = {'err': repr(e)[:200]}; json.dump(res, open(OUT, 'w')); return res[k]
    res[k] = {'score': sc, 'nz': nz, 'ok': ok, 't': round(time.time() - t, 1), 'tag': tag}
    json.dump(res, open(OUT, 'w'))
    print('%-22s score %d  nz %d  ok %s  %.0fs' % (k[:22], sc, nz, ok, time.time() - t), flush=True)
    if sc > 39026:
        json.dump({'x_%d' % i: v[i] for i in range(NV) if v[i]}, open('runs/BEAT_%d_%s.json' % (sc, k[:40]), 'w'))
        print('*** BEATS BASELINE ***', flush=True)
    return res[k]

if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'single'
    do([24601, 2081], 'deliverable-control')
    if mode == 'single':
        for b in sorted(pins, key=int):
            do([int(b)], 'single')
    elif mode == 'pairs':
        rnd = random.Random(7)
        allb = sorted(int(b) for b in pins)
        for _ in range(4000):
            do(rnd.sample(allb, 2), 'pair')
    elif mode == 'triples':
        rnd = random.Random(11)
        allb = sorted(int(b) for b in pins)
        for _ in range(2000):
            do(rnd.sample(allb, 3), 'triple')
