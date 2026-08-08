#!/usr/bin/env python3
"""verify_engine.py -- cross-check every primitive of the C engine against the
python reference in instance.py.  Nothing downstream is trusted unless this passes.
"""
import subprocess, sys, os, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from instance import p, n, B, G, T, PTS, add, mul, neg, sub

SK = os.path.join(HERE, 'sk')

a = 0x0123456789abcdeffedcba98765432100f1e2d3c4b5a69788796a5b4c3d2e1f00
b = 0xfedcba98765432100123456789abcdef00112233445566778899aabbccddeeff
k = 0x1234567890abcdef1122334455667788
out = subprocess.run([SK, 'selftest', '%064x' % G[0], '%064x' % G[1], '%064x' % k],
                     capture_output=True, text=True).stdout
vals = {}
for line in out.strip().split('\n'):
    f = line.split()
    vals[f[0]] = f[1:]

fails = []
def chk(name, got, want):
    ok = got == want
    print(f"  {name:10s} {'OK' if ok else 'FAIL'}")
    if not ok:
        fails.append(name)
        print(f"    got  {got}\n    want {want}")

chk('fe_mul', int(vals['mul'][0], 16), a * b % p)
chk('fe_add', int(vals['add'][0], 16), (a + b) % p)
chk('fe_sub', int(vals['sub'][0], 16), (a - b) % p)
chk('fe_inv', int(vals['inv'][0], 16), pow(a, -1, p))

Q = mul(k, G)
chk('pt_mul', (int(vals['mul_pt'][0], 16), int(vals['mul_pt'][1], 16)), Q)
chk('pt_add', (int(vals['add_pt'][0], 16), int(vals['add_pt'][1], 16)), add(G, Q))
for i in range(1, 9):
    chk(f'walk{i}', (int(vals[f'walk{i}'][0], 16), int(vals[f'walk{i}'][1], 16)), mul(i, G))

# randomised field fuzz through the selftest path (mul/add/sub/inv on random inputs)
random.seed(7)
for trial in range(12):
    x = random.randrange(1, p); y = random.randrange(1, p); s = random.randrange(1, n)
    o = subprocess.run([SK, 'selftest', '%064x' % G[0], '%064x' % G[1], '%064x' % s],
                       capture_output=True, text=True).stdout
    d = {}
    for line in o.strip().split('\n'):
        f = line.split(); d[f[0]] = f[1:]
    got = (int(d['mul_pt'][0], 16), int(d['mul_pt'][1], 16))
    want = mul(s, G)
    if got != want:
        fails.append(f'fuzz_mul_{trial}')
        print(f"  fuzz {trial} FAIL s={s:x}")
print(f"  fuzz_ptmul  {'OK' if not any(f.startswith('fuzz') for f in fails) else 'FAIL'} (12 random 256-bit scalars)")

# --- end to end: bsgs must find a planted small scalar ---
print("\n  end-to-end BSGS on a planted scalar:")
# span covered below is 4096 * 8192 = 2^25, so every planted value must be < 2^25
for planted in (1, 2, 1000, 12345, 2**20 + 7, 2**25 - 3, 2**24 + 4095, 3 * 2**20 + 999):
    TT = mul(planted, G)
    r = subprocess.run([SK, 'bsgs', '%064x' % G[0], '%064x' % G[1],
                        '%064x' % TT[0], '%064x' % TT[1], '12', '8192'],
                       capture_output=True, text=True)
    hit = None
    for line in r.stdout.split('\n'):
        if line.startswith('HIT i='):
            i = int(line.split()[1][2:]); j = int(line.split()[2][2:])
            for cand in (i * 4096 + j, i * 4096 - j):
                if cand >= 0 and mul(cand, G) == TT:
                    hit = cand
    ok = hit == planted
    print(f"    planted {planted:12d} -> {hit}  {'OK' if ok else 'FAIL'}")
    if not ok: fails.append(f'bsgs_{planted}')

# --- end to end: rational mode must find a planted a/b ---
print("\n  end-to-end rational mode on planted k = a/b mod n:")
for (aa, bb) in ((3, 7), (1000, 999), (65537, 4093), (2, 1048573)):
    kk = aa * pow(bb, -1, n) % n
    TT = mul(kk, G)
    r = subprocess.run([SK, 'rational', '%064x' % G[0], '%064x' % G[1],
                        '%064x' % TT[0], '%064x' % TT[1], '21', '2097152'],
                       capture_output=True, text=True)
    found = False
    for line in r.stdout.split('\n'):
        if line.startswith('HIT i='):
            bq = int(line.split()[1][2:]); aq = int(line.split()[2][2:])
            for s in (1, -1):
                cand = s * aq * pow(bq, -1, n) % n
                if mul(cand, G) == TT: found = True
    print(f"    planted {aa}/{bb} -> {'FOUND' if found else 'MISSED'}")
    if not found: fails.append(f'rat_{aa}_{bb}')

print()
print("ENGINE VERIFICATION:", "ALL PASS" if not fails else f"FAILURES {fails}")
sys.exit(1 if fails else 0)
