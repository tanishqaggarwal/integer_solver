#!/usr/bin/env python3
"""Verify, for ALL 96 stages, that the stage's three checks are satisfied exactly when the stage's
output pair equals chordK(inputA, inputB) with the one universal constant K.

The three checked wires of a stage depend only on that stage's six free inputs, so the law can be probed
directly by setting those six inputs -- no gate activation and no repair needed.
Roles (which two of the three pairs are inputs, which is the output) are discovered by testing in which
pair the checks are LINEAR with an invertible 2x2 Jacobian.
"""
import sys, os, json, itertools, random, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fwd import Engine, NV
E = Engine()
p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
K = 97553848499418123410591666447050222001188385549510401465815187079080512838891
tree = json.load(open(os.path.join(HERE, 'tree96.json')))

# map each stage gate -> its three checked wires
import re, collections
gates = collections.defaultdict(list)
for a in E.res:
    m = re.match(r'^\(\(x(\d+)\*x(\d+)\)\+x(\d+)\)$', a)
    if m: gates[int(m.group(1))].append(int(m.group(2))); continue
    m = re.match(r'^\(\((\d+)\*\(x(\d+)\*x(\d+)\)\)-x(\d+)\)$', a)
    if m: gates[int(m.group(2))].append(int(m.group(3))); continue
    m = re.match(r'^\(\(x(\d+)\*x(\d+)\)-\((\d+)\*x(\d+)\)\)$', a)
    if m: gates[int(m.group(1))].append(int(m.group(2))); continue
stages = {g: v for g, v in gates.items() if len(v) >= 3}
resi = {a: i for i, a in enumerate(E.res)}
wire_atom = {}
for g, ws in stages.items():
    wire_atom[g] = ws


def chordK(A, B):
    ax, ay = A; bx, by = B
    if (bx - ax) % p == 0: return None
    l = (by - ay) * pow(bx - ax, p - 2, p) % p
    x = (l * l - ax - bx - K) % p
    return (x, (l * (ax - x) - ay) % p)


def probe(six, ws, assign):
    v = [0] * NV
    for k, x in assign.items(): v[k] = x
    E.run(v)
    return [v[w] % p for w in ws]


def try_roles(six, ws, rnd):
    """Search ALL role splits; return a chordK match if any exists."""
    fallback = None
    for out in itertools.combinations(range(6), 2):
        rest = [i for i in range(6) if i not in out]
        for inA in itertools.combinations(rest, 2):
            inB = tuple(i for i in rest if i not in inA)
            A = (rnd.randrange(1, p), rnd.randrange(1, p))
            B = (rnd.randrange(1, p), rnd.randrange(1, p))
            base = {}
            for idx, val in zip(inA, A): base[six[idx]] = val
            for idx, val in zip(inB, B): base[six[idx]] = val
            for idx in out: base[six[idx]] = 0
            b0 = probe(six, ws, base)
            c1 = [(x - y) % p for x, y in zip(probe(six, ws, {**base, six[out[0]]: 1}), b0)]
            c2 = [(x - y) % p for x, y in zip(probe(six, ws, {**base, six[out[1]]: 1}), b0)]
            det = (c1[0] * c2[1] - c1[1] * c2[0]) % p
            if det == 0: continue
            di = pow(det, p - 2, p)
            dx = ((-b0[0]) * c2[1] + b0[1] * c2[0]) % p * di % p
            dy = (c1[0] * (-b0[1]) + c1[1] * b0[0]) % p * di % p
            if (c1[2] * dx + c2[2] * dy + b0[2]) % p: continue      # third check inconsistent
            got = (dx, dy)
            for P, Q in ((A, B), (B, A)):
                if chordK(P, Q) == got:
                    return (inA, inB, out, 'chordK', (P is A))
            if fallback is None:
                fallback = (inA, inB, out, 'LINEAR_BUT_NOT_chordK', got)
    return fallback


if __name__ == '__main__':
    rnd = random.Random(5)
    ok = 0; other = []; fail = []; skipped = []
    t0 = time.time()
    for g in sorted(stages, key=lambda g: -len(tree[str(g)]['gsup'])):
        six = tree[str(g)]['six']; ws = wire_atom[g]
        if len(six) != 6 or len(ws) < 3:
            skipped.append((g, len(six), len(ws))); continue
        r = try_roles(six, ws, rnd)
        if r is None: fail.append(g); continue
        if r[3] == 'chordK': ok += 1
        else: other.append((g, r[3]))
    print('stages verified to obey chordK exactly : %d of %d' % (ok, len(stages)))
    print('stages linear but NOT chordK           : %d  %s' % (len(other), other[:6]))
    print('stages with no consistent role split   : %d  %s' % (len(fail), fail[:12]))
    print('stages skipped (support != 6 free inputs) : %d  %s' % (len(skipped), skipped[:8]))
    print('elapsed %.0fs' % (time.time() - t0))
