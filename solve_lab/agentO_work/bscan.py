"""STEP 6 — the scan, run as a MEASUREMENT of structure rather than as a search.

The rate says sampling cannot hit the target.  What a scan can still tell us is whether the
four boundary quantities a configuration would have to move are structured:

   K1 = x_7068 - x_2099      (const of a23616)
   L  = x_4432 - x_19964     (const of a23618)
   K2 = 5113045 * x_9118     (external part of a36660)
   J  = x_8731               (external part of a36662, x_7075 = 1)

Questions: do they move with the configuration at all?  Do they move independently?  How many
distinct residues mod p do they take?  Does the small modulus behave differently from p?
"""
import sys, json, time, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentO_work')
import simO, engine as E, harness as H

OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/bscan.log', 'w', buffering=1)
P = 115792089237316195423570985008687907853269984665640564039457584007908834671663
SMALL = 2458959


def say(*a):
    print(*a, file=LOG)


VARS = [7068, 2099, 4432, 19964, 9118, 8731, 7075]


def quantities(v):
    return {'K1': v[7068] - v[2099], 'L': v[4432] - v[19964],
            'K2': 5113045 * v[9118], 'J': v[7075] * v[8731]}


FE = json.load(open(OD + '/feasbits.json'))
A, B = FE['A'], FE['B']
say('feasible a-bits %d, b-bits %d' % (len(A), len(B)))

base = dict(simO.C.base)
v0 = E.forward(base)
q0 = quantities(v0)
say('\nat E cfg0:')
for k, x in q0.items():
    say('   %-3s = %d bits   mod p = %s   mod %d = %d'
        % (k, abs(x).bit_length(), 'ZERO' if x % P == 0 else '%d bits' % (x % P).bit_length(),
           SMALL, x % SMALL))

rows = []
t0 = time.time()
cfgs = [('empty', [])]
cfgs += [('a%d' % b, [b]) for b in A[:12]]
cfgs += [('b%d' % b, [b]) for b in B[:12]]
cfgs += [('a%d+b%d' % (A[i], B[i]), [A[i], B[i]]) for i in range(10)]
for tag, bits in cfgs:
    s = dict(base)
    for b in bits:
        s[b] = 1
    try:
        v = E.forward(s)
    except Exception as ex:
        say('  %-14s ERR %s' % (tag, type(ex).__name__))
        continue
    q = quantities(v)
    rows.append((tag, q))
    say('  %-14s K1=%4db L=%4db K2=%5db J=%5db | modp zero: K1=%s L=%s K2=%s J=%s'
        % (tag, abs(q['K1']).bit_length(), abs(q['L']).bit_length(),
           abs(q['K2']).bit_length(), abs(q['J']).bit_length(),
           q['K1'] % P == 0, q['L'] % P == 0, q['K2'] % P == 0, q['J'] % P == 0))
say('(%.0fs)' % (time.time() - t0))

say('\n--- variability across %d configurations' % len(rows))
for k in ['K1', 'L', 'K2', 'J']:
    vals = [q[k] for _, q in rows]
    resp = [x % P for x in vals]
    ress = [x % SMALL for x in vals]
    say('  %-3s distinct values %3d/%d   distinct mod p %3d   distinct mod %d %3d   always zero mod p: %s'
        % (k, len(set(vals)), len(vals), len(set(resp)), SMALL, len(set(ress)),
           all(x % P == 0 for x in vals)))

say('\n--- are the four correlated?  (do any two always move together)')
for i, k1 in enumerate(['K1', 'L', 'K2', 'J']):
    for k2 in ['K1', 'L', 'K2', 'J'][i + 1:]:
        pairs = {(q[k1], q[k2]) for _, q in rows}
        d1 = len({q[k1] for _, q in rows})
        d2 = len({q[k2] for _, q in rows})
        say('  %s vs %s: %d distinct pairs (independent would be up to %d)'
            % (k1, k2, len(pairs), d1 * d2))
json.dump([[t, {k: str(x) for k, x in q.items()}] for t, q in rows],
          open(OD + '/bscan.json', 'w'))
say('DONE')
