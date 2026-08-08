#!/usr/bin/env python3
"""h4_algebraic.py -- hypothesis 4 and friends: k as a closed-form / algebraic /
textual constant.  Every candidate is tested by the exact predicate mul(k,G)==T.

Because mul() is the only oracle, we can afford ~1e5 candidates in python
(~4 ms each).  We therefore enumerate deliberately rather than exhaustively and
record the exact families covered.
"""
import sys, os, json, time, itertools, hashlib
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from instance import p, n, B, G, T, PTS, add, mul, neg

t_start = time.time()
TESTED = 0
FOUND = []

def test(k, why):
    """exact check.  returns True on success."""
    global TESTED
    k %= n
    TESTED += 1
    if mul(k, G) == T:
        FOUND.append((k, why))
        print(f"\n*** SOLVED ***  k = {k}   ({why})\n")
        return True
    return False

# ---------------------------------------------------------------- families
def fam_powers():
    """k = c^e mod n for small c,e ; and k = c^e (no reduction) when < 2^256"""
    cnt = 0
    for c in range(2, 512):
        for e in range(1, 300):
            if test(pow(c, e, n), f"{c}^{e} mod n"): return cnt, True
            cnt += 1
    return cnt, False

def fam_shift():
    """k = s * 2^j for small s, any j"""
    cnt = 0
    for s in list(range(1, 4096)):
        for j in range(0, 256):
            v = s << j
            if v >= 1 << 256: break
            if test(v, f"{s}*2^{j}"): return cnt, True
            cnt += 1
    return cnt, False

def fam_ndiv():
    """k = floor(n/c), floor(n*a/b), n +- small, n/2 +- small etc."""
    cnt = 0
    for c in range(1, 20000):
        if test(n // c, f"floor(n/{c})"): return cnt, True
        cnt += 1
    for a in range(1, 200):
        for b in range(a + 1, 200):
            if test(n * a // b, f"floor(n*{a}/{b})"): return cnt, True
            cnt += 1
    for c in range(1, 20000):
        if test(p // c, f"floor(p/{c})"): return cnt, True
        cnt += 1
    return cnt, False

def fam_repunit():
    """repunits / repdigits in bases 2..16, and all-ones bit masks"""
    cnt = 0
    for base in range(2, 17):
        for d in range(1, base):
            v = 0
            for L in range(1, 400):
                v = v * base + d
                if v >= 1 << 256: break
                if test(v, f"repdigit {d} x{L} base {base}"): return cnt, True
                cnt += 1
    # bit masks: (2^a - 1) << b, and (2^a-1)*2^b + small
    for a in range(1, 257):
        for b in range(0, 257 - a):
            if test(((1 << a) - 1) << b, f"(2^{a}-1)<<{b}"): return cnt, True
            cnt += 1
    return cnt, False

def fam_pow2pm():
    """2^i +- small, and 2^i +- 2^j"""
    cnt = 0
    for i in range(0, 257):
        for d in range(-4096, 4097):
            v = (1 << i) + d
            if v <= 0: continue
            if test(v, f"2^{i}{d:+d}"): return cnt, True
            cnt += 1
    return cnt, False

def fam_ascii():
    """k = big-endian / little-endian bytes of short ASCII strings"""
    words = []
    base = ["password","secret","satoshi","bitcoin","hello","hello world","test","key",
            "private","privatekey","anthropic","claude","annealer","anneal","qubo",
            "secp256k1","ecdlp","dlog","flag","admin","root","abc","abcdef","1234",
            "12345678","deadbeef","cafebabe","0","1","a","aaa","the answer","42",
            "correct horse battery staple","letmein","opensesame","xyzzy","foobar",
            "integer_solver","solve_lab","EQUATIONS","equations.txt","puzzle",
            "quantum","annealing","simulated annealing","discrete log","structured",
            "low entropy","easy","gotcha","not so fast","nice try","brute force"]
    for w in base:
        words.append(w); words.append(w.upper()); words.append(w.capitalize())
    cnt = 0
    seen = set()
    for w in words:
        if w in seen: continue
        seen.add(w)
        bs = w.encode()
        if len(bs) > 32: bs = bs[:32]
        for order in ('big', 'little'):
            if test(int.from_bytes(bs, order), f"ascii {w!r} {order}"): return cnt, True
            cnt += 1
        # sha256 of the word -- the single most common way a demo scalar is made
        h = hashlib.sha256(bs).digest()
        for order in ('big', 'little'):
            if test(int.from_bytes(h, order), f"sha256({w!r}) {order}"): return cnt, True
            cnt += 1
        h = hashlib.sha256(w.encode()).hexdigest()
        cnt += 0
    return cnt, False

def fam_decimal_digits():
    """digits of pi, e, sqrt2, phi, ln2 as an integer, various lengths and bases"""
    from sympy import pi as spi, E as sE, sqrt, Rational, N, GoldenRatio, log
    import sympy
    cnt = 0
    consts = {'pi': spi, 'e': sE, 'sqrt2': sqrt(2), 'phi': GoldenRatio,
              'ln2': log(2), 'sqrt3': sqrt(3), 'sqrt5': sqrt(5), 'gamma': sympy.EulerGamma,
              'catalan': sympy.Catalan}
    for name, c in consts.items():
        s = sympy.N(c, 120)
        digits = str(s).replace('.', '').replace('-', '')
        for L in range(1, 78):
            v = int(digits[:L])
            if v >= 1 << 256: break
            if test(v, f"first {L} decimal digits of {name}"): return cnt, True
            cnt += 1
        # fractional part scaled to 2^256 (the "nothing up my sleeve" construction)
        frac = sympy.N(c - sympy.floor(c), 100)
        for bits in (128, 160, 192, 224, 251, 255, 256):
            v = int(sympy.floor(frac * 2**bits))
            if test(v, f"frac({name})*2^{bits}"): return cnt, True
            cnt += 1
    return cnt, False

def fam_palindrome_small():
    """binary palindromes of low complexity: 2^a+2^(255-a) style, and decimal palindromes"""
    cnt = 0
    for a in range(0, 256):
        for L in (256, 255, 254, 128, 64, 32):
            if a >= L: continue
            v = (1 << a) | (1 << (L - 1 - a))
            if test(v, f"palindromic pair 2^{a}+2^{L-1-a}"): return cnt, True
            cnt += 1
    return cnt, False

def fam_factorial_fib():
    cnt = 0
    f = 1
    for i in range(1, 200):
        f *= i
        if f < 1 << 256:
            if test(f, f"{i}!"): return cnt, True
            cnt += 1
        if test(f % n, f"{i}! mod n"): return cnt, True
        cnt += 1
    a, b = 0, 1
    for i in range(1, 400):
        a, b = b, a + b
        if test(a % n, f"fib({i})"): return cnt, True
        cnt += 1
    # primorials, mersennes, fermats
    from sympy import prime, primorial, isprime
    for i in range(1, 60):
        if test(int(primorial(i)) % n, f"primorial({i})"): return cnt, True
        cnt += 1
    for i in range(1, 300):
        if test((1 << i) - 1, f"2^{i}-1"): return cnt, True
        cnt += 1
    return cnt, False

def fam_curve_constants():
    """k built out of the instance's own constants"""
    cnt = 0
    cands = {}
    for name, v in (('p', p), ('n', n), ('B', B), ('Gx', G[0]), ('Gy', G[1]),
                    ('Tx', T[0]), ('Ty', T[1])):
        cands[name] = v
    base = dict(cands)
    for name, v in base.items():
        for d in range(-2048, 2049):
            if test(v + d, f"{name}{d:+d}"): return cnt, True
            cnt += 1
        for c in range(2, 64):
            if test(v // c, f"{name}//{c}"): return cnt, True
            if test(v * c % n, f"{name}*{c}"): return cnt, True
            if test(pow(v, c, n), f"{name}^{c} mod n"): return cnt, True
            cnt += 3
        if test(pow(v, -1, n), f"{name}^-1 mod n"): return cnt, True
        cnt += 1
    for a, b in itertools.permutations(base, 2):
        for op, f in (('+', lambda x,y: x+y), ('-', lambda x,y: x-y),
                      ('*', lambda x,y: x*y), ('^', lambda x,y: x ^ y)):
            if test(f(base[a], base[b]) % n, f"{a}{op}{b}"): return cnt, True
            cnt += 1
        if test(base[a] * pow(base[b], -1, n) % n, f"{a}/{b}"): return cnt, True
        cnt += 1
    return cnt, False

def fam_hash_of_points():
    """k = H(G) style derivations"""
    cnt = 0
    blobs = {
        'Gx': G[0].to_bytes(32,'big'), 'Gy': G[1].to_bytes(32,'big'),
        'Tx': T[0].to_bytes(32,'big'), 'Ty': T[1].to_bytes(32,'big'),
        'G':  G[0].to_bytes(32,'big') + G[1].to_bytes(32,'big'),
        'T':  T[0].to_bytes(32,'big') + T[1].to_bytes(32,'big'),
        'B':  B.to_bytes(32,'big'), 'p': p.to_bytes(32,'big'), 'n': n.to_bytes(32,'big'),
        'Gc': bytes([2 + (G[1] & 1)]) + G[0].to_bytes(32,'big'),
        'Tc': bytes([2 + (T[1] & 1)]) + T[0].to_bytes(32,'big'),
    }
    algos = ['sha256','sha512','sha1','md5','sha3_256','blake2b','blake2s','sha224','sha384']
    for bn, bv in blobs.items():
        for al in algos:
            try: h = hashlib.new(al, bv).digest()
            except Exception: continue
            for order in ('big','little'):
                if test(int.from_bytes(h[:32], order), f"{al}({bn}) {order}"): return cnt, True
                cnt += 1
            # double hash
            h2 = hashlib.new(al, h).digest()
            if test(int.from_bytes(h2[:32], 'big'), f"{al}^2({bn})"): return cnt, True
            cnt += 1
    return cnt, False

FAMILIES = [
    ('c^e mod n, 2<=c<512, 1<=e<300', fam_powers),
    ('s*2^j, 1<=s<4096, 0<=j<256', fam_shift),
    ('floor(n/c) and floor(p/c), c<20000; floor(n*a/b), a<b<200', fam_ndiv),
    ('repdigits base 2..16 and (2^a-1)<<b for all a,b', fam_repunit),
    ('2^i + d, 0<=i<=256, |d|<=4096', fam_pow2pm),
    ('ASCII words / sha256 of words (~170 strings, both endiannesses)', fam_ascii),
    ('decimal digit prefixes and 2^b-scaled fractions of pi,e,phi,sqrt2,ln2,gamma,Catalan', fam_decimal_digits),
    ('binary palindromic pairs 2^a+2^(L-1-a)', fam_palindrome_small),
    ('factorials, fibonacci, primorials, 2^i-1', fam_factorial_fib),
    ('instance constants p,n,B,Gx,Gy,Tx,Ty +-2048, /c, *c, ^c, pairwise ops', fam_curve_constants),
    ('cryptographic hashes of the instance points/constants', fam_hash_of_points),
]

if __name__ == '__main__':
    results = []
    for label, fn in FAMILIES:
        t0 = time.time()
        cnt, hit = fn()
        dt = time.time() - t0
        results.append({'family': label, 'candidates': cnt, 'seconds': round(dt, 1), 'hit': hit})
        print(f"{'HIT ' if hit else 'none'}  {cnt:9d} cand  {dt:7.1f}s   {label}")
        if hit: break
    total = time.time() - t_start
    print(f"\ntotal exact mul(k,G)==T evaluations: {TESTED}   wall {total:.1f}s")
    json.dump({'results': results, 'total_tested': TESTED, 'found': [[str(k), w] for k, w in FOUND],
               'wall_seconds': round(total, 1)},
              open(os.path.join(HERE, 'h4_algebraic.json'), 'w'), indent=1)
    print("FOUND:", FOUND if FOUND else "nothing")
