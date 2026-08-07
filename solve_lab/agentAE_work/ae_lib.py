#!/usr/bin/env python3
"""agent AE -- shared curve library + kangaroo driver."""
import json, os, random, subprocess, math, sys

HERE = os.path.dirname(os.path.abspath(__file__))
_d = json.load(open(os.path.join(HERE, 'ae_data.json')))
p = int(_d['p']); a = int(_d['a']); b = int(_d['b']); N = int(_d['N'])
G = (int(_d['G'][0]), int(_d['G'][1]))
T = (int(_d['T'][0]), int(_d['T'][1]))
beta = int(_d['beta']); lam = int(_d['lam'])
ladder = [(int(x), int(y)) for x, y in _d['ladder']]

def inv(z): return pow(z, p - 2, p)
def add(P, Q):
    if P is None: return Q
    if Q is None: return P
    if P[0] == Q[0]:
        if (P[1] + Q[1]) % p == 0: return None
        l = (3 * P[0] * P[0]) * inv(2 * P[1]) % p
    else:
        l = (Q[1] - P[1]) * inv(Q[0] - P[0]) % p
    x = (l * l - P[0] - Q[0]) % p
    return (x, (l * (P[0] - x) - P[1]) % p)
def neg(P): return None if P is None else (P[0], (-P[1]) % p)
def sub(P, Q): return add(P, neg(Q))
def mul(k, P):
    k %= N
    R = None; Q = P
    while k:
        if k & 1: R = add(R, Q)
        Q = add(Q, Q); k >>= 1
    return R
def mulG(k):
    k %= N; R = None; i = 0
    while k:
        if k & 1: R = add(R, ladder[i])
        k >>= 1; i += 1
    return R
def phi(P):  # = [lam]P
    return None if P is None else ((beta * P[0]) % p, P[1])
def oncurve(P):
    return P is None or (P[1] * P[1] - P[0] ** 3 - b) % p == 0

def hx(v): return format(v, 'x')

def gen_input(path, Q, R, seed, nj=32):
    """Write an engine input file solving log_G(Q) in [0,2^R)."""
    rng = random.Random(seed)
    m = 1 << (R // 2)                      # mean jump ~ sqrt(L)
    lines = []
    lines.append('%s %s' % (hx(G[0]), hx(G[1])))
    lines.append('%s %s' % (hx(Q[0]), hx(Q[1])))
    lines.append(str(R))
    lines.append(str(nj))
    jds = []
    for i in range(nj):
        jd = 1 + rng.randrange(2 * m)
        J = mulG(jd)
        assert J is not None
        jds.append(jd)
        lines.append('%s %s %s' % (hx(jd), hx(J[0]), hx(J[1])))
    for L in ladder:
        lines.append('%s %s' % (hx(L[0]), hx(L[1])))
    open(path, 'w').write('\n'.join(lines) + '\n')
    return sum(jds) / float(nj)

def run_kangaroo(tag, Q, R, threads=2, kpt=1024, dpbits=None, log2max=None,
                 seed=1, tablebits=21, workdir=None, quiet=False):
    """Returns dict with outcome.  Q must satisfy log_G(Q) in [0,2^R) for a hit."""
    wd = workdir or HERE
    inp = os.path.join(wd, 'in_%s.txt' % tag)
    meanj = gen_input(inp, Q, R, seed)
    K = threads * kpt
    if dpbits is None:
        # per-kangaroo jumps ~ 2*sqrt(L)/K ; want >= 64 DP intervals each
        perk = max(1.0, 2.0 * (2 ** (R / 2.0)) / K)
        dpbits = max(0, int(math.log2(perk / 64.0)))
    if log2max is None:
        log2max = int(R / 2.0 + 3)         # 8x the ~2*sqrt(L) expectation
    cmd = [os.path.join(HERE, 'aekang'), inp, str(threads), str(kpt),
           str(dpbits), str(log2max), str(seed), str(tablebits)]
    pr = subprocess.run(cmd, capture_output=True, text=True)
    out = pr.stdout; err = pr.stderr
    res = dict(tag=tag, R=R, rc=pr.returncode, meanjump=meanj, dpbits=dpbits,
               log2max=log2max, K=K, cands=[], done=None, stderr_tail=err[-400:])
    for line in out.splitlines():
        if line.startswith('CAND'):
            f = line.split()
            st = int(f[1][3:]); sw = int(f[2][3:])
            res['cands'].append((st, sw))
        elif line.startswith('DONE'):
            res['done'] = line
            for kv in line.split()[1:]:
                k, v = kv.split('=')
                try: res[k] = float(v) if '.' in v else int(v)
                except ValueError: res[k] = v
    if not quiet:
        print('[%s] rc=%d %s' % (tag, pr.returncode, res.get('done')))
        if pr.returncode != 0: print('  stderr:', err[-300:])
    return res

def check_cands(res, Q):
    """Turn engine candidates into verified discrete logs of Q base G."""
    good = []
    for st, sw in res['cands']:
        for cand in ((st - sw) % N, (-st - sw) % N):
            if mulG(cand) == Q:
                good.append(cand)
    return good
