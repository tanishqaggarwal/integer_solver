#!/usr/bin/env python3
"""agent AE -- driver + validator for the small-quotient sweep aequot.

Family:  k0 == +/- a * lam^-e * b^-1 (mod N),  1 <= a,b <= 2^LOGM, e in {0,1,2}.
Contains as special cases b=1 (k0 small), a=1 (k0 = 1/b), and every "nice
rational" scalar.  A miss excludes the whole family exhaustively -- unlike the
kangaroo this instrument is deterministic.
"""
import sys, os, json, subprocess, time
import ae_lib as L

N = L.N

def gen(path, Tpt, LOGM, W):
    M = 1 << LOGM
    assert M % W == 0
    chunk = M // W
    lines = ['%s %s' % (L.hx(L.G[0]), L.hx(L.G[1])),
             '%s %s' % (L.hx(Tpt[0]), L.hx(Tpt[1])),
             L.hx(L.beta), '%d %d' % (LOGM, W)]
    for base in (L.G, Tpt):
        D = L.mul(chunk, base)
        P = base           # j=0 -> scalar 1
        pts = []
        for j in range(W):
            pts.append(P)
            P = L.add(P, D)
        # verify two of them the slow way
        for j in (1, W - 1):
            assert pts[j] == L.mul(j * chunk + 1, base), 'lane start check failed'
        for P2 in pts:
            lines.append('%s %s' % (L.hx(P2[0]), L.hx(P2[1])))
    open(path, 'w').write('\n'.join(lines) + '\n')
    return M

def run(tag, Tpt, LOGM, W=1024):
    inp = 'q_in_%s.txt' % tag; outp = 'q_out_%s.txt' % tag
    M = gen(inp, Tpt, LOGM, W)
    t0 = time.time()
    pr = subprocess.run([os.path.join(L.HERE, 'aequot'), inp, outp],
                        capture_output=True, text=True)
    el = time.time() - t0
    matches = []; done = None
    for line in open(outp):
        if line.startswith('MATCH'):
            f = dict(kv.split('=') for kv in line.split()[1:])
            matches.append((int(f['a']), int(f['b']), int(f['e'])))
        elif line.startswith('DONE'):
            done = line.strip()
    os.remove(inp)
    return dict(tag=tag, M=M, rc=pr.returncode, secs=round(el, 1), done=done,
                matches=matches, stderr=pr.stderr.strip())

def verify(matches, Tpt):
    """turn (a,b,e) into verified scalars k with k*G == Tpt"""
    good = []
    for (a, b, e) in matches:
        for s in (1, -1):
            k = (s * a * pow(L.lam, -e, N) * pow(b, -1, N)) % N
            if L.mulG(k) == Tpt:
                good.append((a, b, e, s, k))
    return good

if __name__ == '__main__':
    mode = sys.argv[1]
    LOGM = int(sys.argv[2]) if len(sys.argv) > 2 else 22
    if mode == 'plant':
        # three plants: plain quotient, lam-twisted, negated
        rng = __import__('random').Random(4242)
        allok = True
        for i, (e0, s0) in enumerate([(0, 1), (1, 1), (2, -1), (0, -1)]):
            a0 = rng.randrange(1, 1 << LOGM); b0 = rng.randrange(1, 1 << LOGM)
            k = (s0 * a0 * pow(L.lam, -e0, N) * pow(b0, -1, N)) % N
            Tp = L.mulG(k)
            r = run('plant%d' % i, Tp, LOGM)
            g = verify(r['matches'], Tp)
            found = any(a == a0 and b == b0 for (a, b, e, s, kk) in g)
            recovered = any(kk == k for (a, b, e, s, kk) in g)
            print('plant %d (e=%d,s=%+d) a=%d b=%d : %d matches, exact pair found=%s, k recovered=%s  %s'
                  % (i, e0, s0, a0, b0, len(r['matches']), found, recovered, r['done']))
            allok &= (found and recovered)
        print('PLANT VALIDATION:', 'PASS' if allok else 'FAIL')
    elif mode == 'real':
        r = run('real', L.T, LOGM)
        g = verify(r['matches'], L.T)
        print(r['done'])
        print('rc=%d  raw matches=%d  verified scalars=%d' % (r['rc'], len(r['matches']), len(g)))
        for x in g: print('  VERIFIED', x)
        json.dump(dict(done=r['done'], rc=r['rc'], nmatch=len(r['matches']),
                       verified=[[a, b, e, s, str(k)] for (a, b, e, s, k) in g]),
                  open('res_quot_%d.json' % LOGM, 'w'), indent=1)
