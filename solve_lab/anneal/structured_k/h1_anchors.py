#!/usr/bin/env python3
"""h1_anchors.py -- hypothesis 1b: k sits within a +-2^SPANBITS/2 window of a
"round" anchor value (n, n/2, 2^255, 2^256-2^32, ...).

For anchor A we hand the C engine the target  T - (A - span/2)*G , so a small
dlog there means k = A - span/2 + (that dlog).  One baby table serves every
anchor.  A planted control anchor is included so a silent failure is impossible.
"""
import sys, os, subprocess, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from instance import p, n, B, G, T, PTS, add, mul, neg, sub

SK = os.path.join(HERE, 'sk')
BB = 24                       # baby table 2^24
GC = 1 << 21                  # giant steps
SPAN = (1 << BB) * GC         # 2^45
HALF = SPAN // 2

anchors = {}
def A(name, v): anchors[name] = v % n

A('0', 0)
A('n', n); A('n-1', n - 1); A('n+1', n + 1)
A('n/2', n // 2); A('(n+1)/2', (n + 1) // 2); A('n/3', n // 3); A('2n/3', 2 * n // 3)
A('n/4', n // 4); A('3n/4', 3 * n // 4)
A('p', p); A('p/2', p // 2)
A('2^255', 1 << 255); A('2^254', 1 << 254); A('2^253', 1 << 253)
A('2^256-2^32', (1 << 256) - (1 << 32)); A('2^256-1', (1 << 256) - 1)
A('2^256-2^32-977', (1 << 256) - (1 << 32) - 977)
for b in (32, 48, 64, 96, 128, 160, 192, 224, 240, 248, 250, 251, 252, 255, 256):
    A(f'2^{b}', 1 << b)
    A(f'2^{b}-1', (1 << b) - 1)
A('B', B); A('Gx', G[0]); A('Gy', G[1]); A('Tx', T[0]); A('Ty', T[1])
A('sqrt(n)', int(n ** 0.5))
A('n-2^32', n - (1 << 32)); A('n-2^45', n - (1 << 45))
# every 1/16th of the group -- a coarse net over the whole keyspace
for i in range(1, 16):
    A(f'{i}n/16', i * n // 16)

# planted control: a value we KNOW is inside a window, to prove the pipeline fires
CONTROL_ANCHOR = (1 << 255)
CONTROL_K = CONTROL_ANCHOR + 123456789
anchors['CONTROL'] = CONTROL_ANCHOR

def build_input():
    lines = []
    meta = []
    for name, a in anchors.items():
        base = (a - HALF) % n
        if name == 'CONTROL':
            tgt = sub(mul(CONTROL_K, G), mul(base, G))
        else:
            tgt = sub(T, mul(base, G))
        if tgt is None:
            print(f"  anchor {name}: T == {a}*G exactly?!")
            if mul(a, G) == T:
                print(f"*** SOLVED *** k = {a}")
                json.dump({'solved': str(a)}, open(os.path.join(HERE, 'SOLVED.json'), 'w'))
                sys.exit(0)
            continue
        lines.append(f"{name} {tgt[0]:064x} {tgt[1]:064x}")
        meta.append((name, a, base))
    return "\n".join(lines) + "\n", meta

if __name__ == '__main__':
    t0 = time.time()
    inp, meta = build_input()
    print(f"{len(meta)} anchors, window +-2^{SPAN.bit_length()-2} around each "
          f"(span 2^{SPAN.bit_length()-1} = {SPAN})")
    basemap = {name: base for name, a, base in meta}
    r = subprocess.run([SK, 'bsgsmulti', f"{G[0]:064x}", f"{G[1]:064x}", str(BB), str(GC)],
                       input=inp, capture_output=True, text=True)
    print(r.stderr.strip())
    hits, results = [], []
    pend = None
    for line in r.stdout.split('\n'):
        f = line.split()
        if not f: continue
        if f[0] in ('MHIT', 'MNOHIT'):
            results.append({'anchor': f[1], 'hit': f[0] == 'MHIT'})
            pend = f[1]
        elif f[0] == 'HITDATA':
            name = f[1]; i = int(f[2][2:]); j = int(f[3][2:])
            M = 1 << BB
            base = basemap[name]
            for s in (1, -1):
                k = (base + i * M + s * j) % n
                if mul(k, G) == T:
                    hits.append((name, k))
                    print(f"\n*** SOLVED ***  k = {k}   (anchor {name})\n")
                elif name == 'CONTROL' and (base + i*M + s*j) % n == CONTROL_K % n:
                    print(f"  control anchor fired correctly (recovered {CONTROL_K})")
    nhit = sum(1 for x in results if x['hit'])
    print(f"anchors probed: {len(results)}, raw x-coord hits: {nhit}, verified solutions: {len(hits)}")
    print(f"wall {time.time()-t0:.1f}s")
    json.dump({'span': SPAN, 'span_bits': SPAN.bit_length() - 1, 'babybits': BB,
               'giantcount': GC, 'anchors': {k: str(v) for k, v in anchors.items()},
               'results': results, 'solutions': [[a, str(k)] for a, k in hits],
               'wall_seconds': round(time.time() - t0, 1)},
              open(os.path.join(HERE, 'h1_anchors.json'), 'w'), indent=1)
    if hits:
        json.dump({'solved': str(hits[0][1])}, open(os.path.join(HERE, 'SOLVED.json'), 'w'))
