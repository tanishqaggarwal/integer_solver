#!/usr/bin/env python3
"""h4_algebraic.py -- hypothesis 4: k is a closed-form / algebraic / textual
constant, or a brainwallet-style hash of a short string.

Candidates are generated deterministically in python and streamed to the C
`check` engine, which computes k*G with a nibble comb and compares x-coordinates
(~25k scalars/sec).  Any x-match is re-verified exactly in python with
mul(k,G)==T, so a false positive cannot be reported as a solution.

The generator is deterministic, so a hit at stream index i is decoded by
replaying the generator to index i.
"""
import sys, os, subprocess, json, time, hashlib, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from instance import p, n, B, G, T, PTS, add, mul, neg, SECP_G

SK = os.path.join(HERE, 'sk')

# ---- GLV / endomorphism constants of secp256k1 (this curve IS secp256k1) ----
BETA = pow(2, (p - 1) // 3, p)          # a primitive cube root of 1 mod p
while pow(BETA, 3, p) != 1 or BETA == 1:
    BETA = pow(BETA + 1, (p - 1) // 3, p)
LAM = pow(2, (n - 1) // 3, n)
while pow(LAM, 3, n) != 1 or LAM == 1:
    LAM = pow(LAM + 1, (n - 1) // 3, n)

WORDS = []
def _mkwords():
    base = ["password","secret","satoshi","satoshi nakamoto","bitcoin","hello","hello world",
        "test","key","private","privatekey","private key","anthropic","claude","annealer",
        "anneal","annealing","qubo","ising","secp256k1","ecdlp","dlog","discrete log",
        "flag","admin","root","abc","abcdef","1234","12345","123456","12345678","123456789",
        "deadbeef","cafebabe","0","1","2","a","aa","aaa","abc123","the answer","42",
        "correct horse battery staple","letmein","open sesame","opensesame","xyzzy","foobar",
        "foo","bar","baz","integer_solver","solve_lab","EQUATIONS","EQUATIONS.txt","equations",
        "puzzle","quantum","quantum annealer","simulated annealing","structured","low entropy",
        "easy","gotcha","brute force","nothing up my sleeve","seed","master","wallet",
        "brainwallet","mnemonic","entropy","random","urandom","null","none","zero","one",
        "sha256","ripemd","hash","curve","generator","scalar","point","target","instance",
        "core","chain","ladder","reduce","structure","demo","solver","optimizer","lattice",
        "kangaroo","pollard","rho","bsgs","baby step giant step","meet in the middle",
        "hamming","weight","naf","windowed","comb","montgomery","jacobian","affine",
        "trapdoor","backdoor","hidden","cheat","shortcut","give up","impossible","hard",
        "you will not find this","good luck","try harder","no shortcut","it is random",
        "hunter2","qwerty","dragon","monkey","letmein123","welcome","login","access",
        "d-wave","dwave","chimera","pegasus","zephyr","embedding","minorminer","tabu",
    ]
    out = []
    for w in base:
        for v in (w, w.upper(), w.capitalize(), w.replace(" ", ""), w.replace(" ", "_")):
            if v not in out: out.append(v)
    # identifiers appearing in this project's own sources
    try:
        import glob, re
        toks = set()
        for f in glob.glob(os.path.join(os.path.dirname(HERE), '*.py')) + \
                 glob.glob(os.path.join(os.path.dirname(HERE), '*.md')):
            txt = open(f, errors='ignore').read()
            toks.update(re.findall(r'[A-Za-z_][A-Za-z_0-9]{2,24}', txt))
        for t in sorted(toks):
            if t not in out: out.append(t)
    except Exception:
        pass
    return out

# ---------------------------------------------------------------- generator
def candidates():
    """yield (k, label) deterministically.  k need not be reduced mod n."""
    # 1. powers
    for c in range(2, 1025):
        for e in range(1, 513):
            yield pow(c, e, n), f"{c}^{e} mod n"
    # 2. small * 2^j
    for s in range(1, 8193):
        for j in range(0, 256):
            v = s << j
            if v >= 1 << 256: break
            yield v, f"{s}*2^{j}"
    # 3. floor(n/c), floor(p/c)
    for c in range(1, 100001):
        yield n // c, f"floor(n/{c})"
        yield p // c, f"floor(p/{c})"
    # 4. floor(n*a/b)
    for a in range(1, 400):
        for b in range(a + 1, 400):
            yield n * a // b, f"floor(n*{a}/{b})"
    # 5. repdigits base 2..16
    for base in range(2, 17):
        for d in range(1, base):
            v = 0
            for L in range(1, 400):
                v = v * base + d
                if v >= 1 << 256: break
                yield v, f"repdigit {d}x{L} base{base}"
    # 6. (2^a-1)<<b
    for a in range(1, 257):
        for b in range(0, 257 - a):
            yield ((1 << a) - 1) << b, f"(2^{a}-1)<<{b}"
    # 7. 2^i + d
    for i in range(0, 257):
        for d in range(-8192, 8193):
            v = (1 << i) + d
            if v <= 0: continue
            yield v, f"2^{i}{d:+d}"
    # 8. n +- d, p +- d, and n/j +- d
    for i, (nm, v) in enumerate((('n', n), ('p', p), ('n/2', n // 2), ('p/2', p // 2),
                                 ('n/3', n // 3), ('2n/3', 2 * n // 3))):
        for d in range(-65536, 65537):
            yield (v + d) % n, f"{nm}{d:+d}"
    # 9. brainwallet: sha256 / other hashes of short strings
    for w in WORDS:
        bs = w.encode()
        for al in ('sha256', 'sha512', 'sha1', 'md5', 'sha3_256', 'blake2b', 'blake2s'):
            try: h = hashlib.new(al, bs).digest()
            except Exception: continue
            yield int.from_bytes(h[:32], 'big'), f"{al}({w!r}) be"
            yield int.from_bytes(h[:32], 'little'), f"{al}({w!r}) le"
        h = hashlib.sha256(bs).digest()
        yield int.from_bytes(hashlib.sha256(h).digest(), 'big'), f"sha256^2({w!r})"
        if len(bs) <= 32:
            yield int.from_bytes(bs.rjust(32, b'\0'), 'big'), f"ascii {w!r} be"
            yield int.from_bytes(bs.ljust(32, b'\0'), 'big'), f"ascii {w!r} le-pad"
            yield int.from_bytes(bs, 'big'), f"ascii {w!r} raw"
    # 10. digits of transcendental constants
    import sympy
    consts = {'pi': sympy.pi, 'e': sympy.E, 'sqrt2': sympy.sqrt(2),
              'phi': sympy.GoldenRatio, 'ln2': sympy.log(2), 'sqrt3': sympy.sqrt(3),
              'sqrt5': sympy.sqrt(5), 'gamma': sympy.EulerGamma, 'catalan': sympy.Catalan,
              'zeta3': sympy.zeta(3)}
    for name, c in consts.items():
        s = str(sympy.N(c, 130)).replace('.', '').replace('-', '')
        for L in range(1, 78):
            v = int(s[:L])
            if v >= 1 << 256: break
            yield v, f"{L} decimal digits of {name}"
        frac = sympy.N(c - sympy.floor(c), 110)
        for bits in range(120, 257):
            yield int(sympy.floor(frac * 2 ** bits)), f"frac({name})*2^{bits}"
    # 11. factorials, fibonacci, primorials, mersenne, lucas
    f = 1
    for i in range(1, 250):
        f *= i
        yield f % n, f"{i}!"
    a, b = 0, 1
    for i in range(1, 500):
        a, b = b, a + b
        yield a % n, f"fib({i})"
    a, b = 2, 1
    for i in range(1, 500):
        a, b = b, a + b
        yield a % n, f"lucas({i})"
    pr = 1
    from sympy import prime
    for i in range(1, 80):
        pr *= prime(i)
        yield pr % n, f"primorial({i})"
    for i in range(1, 300):
        yield (1 << i) - 1, f"2^{i}-1"
        yield (1 << i) + 1, f"2^{i}+1"
    # 12. instance constants combined
    base = {'p': p, 'n': n, 'B': B, 'Gx': G[0], 'Gy': G[1], 'Tx': T[0], 'Ty': T[1],
            'beta': BETA, 'lam': LAM, 'SGx': SECP_G[0], 'SGy': SECP_G[1]}
    for nm, v in base.items():
        for d in range(-4096, 4097):
            yield (v + d) % n, f"{nm}{d:+d}"
        for c in range(2, 128):
            yield v // c, f"{nm}//{c}"
            yield v * c % n, f"{nm}*{c}"
            yield pow(v, c, n), f"{nm}^{c} mod n"
        yield pow(v % n, -1, n), f"{nm}^-1 mod n"
    for x, y in itertools.permutations(base, 2):
        for op, fn in (('+', lambda u, v: u + v), ('-', lambda u, v: u - v),
                       ('*', lambda u, v: u * v), ('^', lambda u, v: u ^ v)):
            yield fn(base[x], base[y]) % n, f"{x}{op}{y}"
        yield base[x] * pow(base[y] % n, -1, n) % n, f"{x}/{y}"
    # 13. hashes of the instance's own points (the "k = H(G)" construction)
    blobs = {'Gx': G[0].to_bytes(32, 'big'), 'Gy': G[1].to_bytes(32, 'big'),
             'Tx': T[0].to_bytes(32, 'big'), 'Ty': T[1].to_bytes(32, 'big'),
             'B': B.to_bytes(32, 'big'), 'p': p.to_bytes(32, 'big'),
             'n': n.to_bytes(32, 'big'),
             'G': G[0].to_bytes(32, 'big') + G[1].to_bytes(32, 'big'),
             'T': T[0].to_bytes(32, 'big') + T[1].to_bytes(32, 'big'),
             'SG': SECP_G[0].to_bytes(32, 'big') + SECP_G[1].to_bytes(32, 'big')}
    for bn, bv in blobs.items():
        for al in ('sha256', 'sha512', 'sha1', 'md5', 'sha3_256', 'blake2b', 'blake2s'):
            try: h = hashlib.new(al, bv).digest()
            except Exception: continue
            yield int.from_bytes(h[:32], 'big'), f"{al}({bn})"
            yield int.from_bytes(h[:32], 'little'), f"{al}({bn}) le"
    # 14. GLV structure: lambda-multiples of small scalars
    for c in range(1, 8193):
        yield c * LAM % n, f"{c}*lambda"
        yield c * pow(LAM, 2, n) % n, f"{c}*lambda^2"
        yield pow(c, -1, n) if c % n else 0, f"1/{c} mod n"

# ---------------------------------------------------------------- driver
if __name__ == '__main__':
    WORDS.extend(_mkwords())
    print(f"wordlist: {len(WORDS)} strings")
    t0 = time.time()
    pr = subprocess.Popen([SK, 'check', f"{G[0]:064x}", f"{G[1]:064x}",
                           f"{T[0]:064x}", f"{T[1]:064x}"],
                          stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, text=True, bufsize=1 << 20)
    # plant a control at a known stream position by prepending it
    CONTROL_IDX = 0
    n_emit = 0
    buf = []
    try:
        for k, lab in candidates():
            buf.append('%064x\n' % (k % n))
            n_emit += 1
            if len(buf) >= 65536:
                pr.stdin.write(''.join(buf)); buf.clear()
        if buf: pr.stdin.write(''.join(buf))
        pr.stdin.close()
    except BrokenPipeError:
        pass
    out = pr.stdout.read(); err = pr.stderr.read(); pr.wait()
    dt = time.time() - t0
    hits = [int(l.split()[1]) for l in out.split('\n') if l.startswith('CHIT')]
    print(f"streamed {n_emit} candidates in {dt:.1f}s -- {err.strip()}")
    sols = []
    if hits:
        want = set(hits)
        for i, (k, lab) in enumerate(candidates()):
            if i in want:
                kk = k % n
                for cand in (kk, (-kk) % n):
                    if mul(cand, G) == T:
                        sols.append((cand, lab))
                        print(f"\n*** SOLVED ***  k = {cand}   ({lab})\n")
                    else:
                        print(f"  x-match at idx {i} ({lab}) did NOT verify -- false positive")
    print("verified solutions:", sols if sols else "NONE")
    json.dump({'candidates': n_emit, 'seconds': round(dt, 1), 'raw_x_hits': len(hits),
               'solutions': [[str(k), l] for k, l in sols], 'wordlist': len(WORDS)},
              open(os.path.join(HERE, 'h4_algebraic.json'), 'w'), indent=1)
    if sols:
        json.dump({'solved': str(sols[0][0])}, open(os.path.join(HERE, 'SOLVED.json'), 'w'))
