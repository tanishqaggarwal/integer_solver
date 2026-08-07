"""Agent AD -- extra curve sizes for the SCALING study (n = 8,10,12,14,16,18,20).

Same requirements as ad_curves.py: short Weierstrass over a prime field, PRIME
order N of exactly `n` bits (so the ladder 2^i*G, i=0..n-1, spans the group the
way the real instance's does), one j=0 (CM by sqrt(-3)) curve and one generic
curve at each size.  Written to ad_curves2.json, merged with ad_curves*.json.
"""
import json, os, sys, random
from sympy import isprime
from ad_curves import qr_table, order_of_curve, find_gen, primes_near
from ad_model import Curve

HERE = os.path.dirname(os.path.abspath(__file__))


def one_size(nbits, seed=7):
    rng = random.Random(seed * 1000 + nbits)
    out = []
    for p in primes_near(nbits, 6, cong3=1):
        qr = qr_table(p)
        got0 = got1 = False
        for _ in range(3000):
            if got0 and got1:
                break
            j0 = not got0
            a = 0 if j0 else rng.randrange(1, p)
            b = rng.randrange(1, p)
            if (4 * a ** 3 + 27 * b * b) % p == 0:
                continue
            N = order_of_curve(p, a, b, qr)
            if not (isprime(N) and N.bit_length() == nbits):
                continue
            G = find_gen(p, a, b, N)
            if G is None:
                continue
            out.append(dict(p=p, a=a, b=b, N=N, j0=bool(j0), G=list(G)))
            if j0:
                got0 = True
            else:
                got1 = True
        if got0 and got1:
            break
    return out


def main():
    sizes = [int(x) for x in (sys.argv[1:] or ['10', '14', '18', '20'])]
    out = {}
    for nb in sizes:
        cur = one_size(nb)
        out[str(nb)] = cur
        for c in cur:
            cv = Curve(c['p'], c['a'], c['b'], c['N'], c['G'])
            assert cv.on(cv.G) and cv.mul(c['N'], cv.G) is None
            print('%2d bits: p=%d a=%d b=%d N=%d j0=%s G=%s' %
                  (nb, c['p'], c['a'], c['b'], c['N'], c['j0'], c['G']))
        sys.stdout.flush()
    with open(os.path.join(HERE, 'ad_curves2.json'), 'w') as f:
        json.dump(out, f, indent=1)
    print('written')


if __name__ == '__main__':
    main()
