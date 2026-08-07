#!/usr/bin/env python3
"""Stage-law verification, v2: for every stage, search all role PARTITIONS and, for each, all
coordinate ORDERINGS (which member of a pair is x and which is y, and which input is first).
The linear solve is done once per partition; orderings are then tested arithmetically, so the whole
sweep costs ~135 forward passes per stage.
Reports, per stage, whether the three checks are satisfied exactly when the output equals
chord-with-offset(inA, inB), and what offset that is.
"""
import sys, os, json, itertools, random, time, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fwd import Engine, NV
E = Engine()
p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
KUNIV = 97553848499418123410591666447050222001188385549510401465815187079080512838891
tree = json.load(open(os.path.join(HERE, 'tree96.json')))
import re
gates = collections.defaultdict(list)
for a in E.res:
    m = re.match(r'^\(\(x(\d+)\*x(\d+)\)\+x(\d+)\)$', a)
    if m: gates[int(m.group(1))].append(int(m.group(2))); continue
    m = re.match(r'^\(\((\d+)\*\(x(\d+)\*x(\d+)\)\)-x(\d+)\)$', a)
    if m: gates[int(m.group(2))].append(int(m.group(3))); continue
    m = re.match(r'^\(\(x(\d+)\*x(\d+)\)-\((\d+)\*x(\d+)\)\)$', a)
    if m: gates[int(m.group(1))].append(int(m.group(2))); continue
stages = {g: v for g, v in gates.items() if len(v) >= 3}


def probe(ws, assign):
    v = [0] * NV
    for k, x in assign.items(): v[k] = x
    E.run(v)
    return [v[w] % p for w in ws]


def offset_of(A, B, out):
    """offset K with out == chordK(A,B); None if the y-law does not hold."""
    if (B[0] - A[0]) % p == 0: return None
    l = (B[1] - A[1]) * pow(B[0] - A[0], p - 2, p) % p
    if (l * (A[0] - out[0]) - A[1]) % p != out[1] % p: return None
    return (l * l - A[0] - B[0] - out[0]) % p


def analyse(g, rnd, trials=2):
    six = tree[str(g)]['six']; ws = stages[g]
    if len(six) != 6: return ('skip', len(six))
    hits = collections.Counter()
    for out in itertools.combinations(range(6), 2):
        rest = [i for i in range(6) if i not in out]
        for inA in itertools.combinations(rest, 2):
            inB = tuple(i for i in rest if i not in inA)
            per_trial = []
            ok = True
            for t in range(trials):
                a0, a1 = rnd.randrange(1, p), rnd.randrange(1, p)
                b0v, b1v = rnd.randrange(1, p), rnd.randrange(1, p)
                base = {six[inA[0]]: a0, six[inA[1]]: a1,
                        six[inB[0]]: b0v, six[inB[1]]: b1v,
                        six[out[0]]: 0, six[out[1]]: 0}
                f0 = probe(ws, base)
                c1 = [(x - y) % p for x, y in zip(probe(ws, {**base, six[out[0]]: 1}), f0)]
                c2 = [(x - y) % p for x, y in zip(probe(ws, {**base, six[out[1]]: 1}), f0)]
                det = (c1[0] * c2[1] - c1[1] * c2[0]) % p
                if det == 0: ok = False; break
                di = pow(det, p - 2, p)
                d0 = ((-f0[0]) * c2[1] + f0[1] * c2[0]) % p * di % p
                d1 = (c1[0] * (-f0[1]) + c1[1] * f0[0]) % p * di % p
                if (c1[2] * d0 + c2[2] * d1 + f0[2]) % p: ok = False; break
                found = set()
                for oswap in (0, 1):
                    O = (d1, d0) if oswap else (d0, d1)
                    for asw in (0, 1):
                        A = (a1, a0) if asw else (a0, a1)
                        for bsw in (0, 1):
                            Bp = (b1v, b0v) if bsw else (b0v, b1v)
                            for P, Q in ((A, Bp), (Bp, A)):
                                k = offset_of(P, Q, O)
                                if k is not None: found.add((oswap, asw, bsw, P is A, k))
                per_trial.append(found)
            if not ok or not per_trial: continue
            common = set.intersection(*per_trial) if len(per_trial) > 1 else per_trial[0]
            for tup in common: hits[tup[4]] += 1
    return ('ok', hits)


if __name__ == '__main__':
    rnd = random.Random(31)
    t0 = time.time()
    good = 0; univ = 0; other = collections.Counter(); none = []; skipped = []
    for g in sorted(stages, key=lambda g: -len(tree[str(g)]['gsup'])):
        kind, r = analyse(g, rnd)
        if kind == 'skip': skipped.append((g, r)); continue
        if not r: none.append(g); continue
        good += 1
        ks = set(r)
        if KUNIV in ks: univ += 1
        else: other[tuple(sorted(ks))[0]] += 1
    print('stages with 6 free inputs analysed      : %d' % (len(stages) - len(skipped)))
    print('  admit a chord-with-offset law         : %d' % good)
    print('  and the offset is the UNIVERSAL K     : %d' % univ)
    print('  offset different from K               : %d  %s' % (sum(other.values()), [str(k)[:30] for k in other][:4]))
    print('  no consistent law found               : %d  %s' % (len(none), none[:10]))
    print('stages skipped (fewer than 6 free ins)  : %d' % len(skipped))
    print('elapsed %.0fs' % (time.time() - t0))
