"""W stage 4: does the 3x2 congruence matrix have rank 2 MOD P at every block?
A rank drop mod P would be a genuine third family (N1,N2 not forced to 0)."""
import sys, os, re, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
from collections import Counter
PVAL = 115792089237316195423570985008687907853269984665640564039457584007908834671663
blocks = json.load(open('w_blocks3.json'))
minmin = None; bad = []; rk = Counter(); allmin = []
gcdP = Counter(); acoef = []
for b in blocks:
    M = [(cg['ring']['c1'], cg['ring']['c2']) for cg in b['congs']]
    mins = [M[i][0]*M[j][1] - M[i][1]*M[j][0] for i in range(3) for j in range(i+1, 3)]
    allmin += mins
    nz = [m for m in mins if m != 0]
    rk['rankQ=2' if nz else 'rankQ<2'] += 1
    if not nz: bad.append(('rankQ<2', b['E']))
    # rank mod P: drops iff ALL minors == 0 mod P
    if all(m % PVAL == 0 for m in mins): bad.append(('rankP<2', b['E']))
    for cg in b['congs']:
        r = cg['ring']
        acoef.append((r['aout'], r['chand']))
        gcdP[math.gcd(r['aout'], PVAL), math.gcd(abs(r['chand']), PVAL)] += 1
print('rank over Q:', dict(rk))
print('blocks with rank<2 (Q or mod P):', bad[:10], 'count', len(bad))
print('max |2x2 minor|:', max(abs(m) for m in allmin), ' (P has %d bits)' % PVAL.bit_length())
print('min |nonzero minor|:', min(abs(m) for m in allmin if m), ' zero minors:', sum(1 for m in allmin if m == 0))
print('gcd(a,P),gcd(c,P) classes:', gcdP.most_common(4))
# structural: the 3x2 matrices distinct across blocks?
sig = Counter(tuple(sorted((cg['ring']['c1'], cg['ring']['c2']) for cg in b['congs'])) for b in blocks)
print('distinct 3x2 matrices:', len(sig), 'of', len(blocks))
# handle multipliers
big_c = [ (b['E'], cg['ring']['aout'], cg['ring']['chand']) for b in blocks for cg in b['congs']
          if abs(cg['ring']['chand']) > 1 ]
print('congruences with |c|>1 on the handle:', len(big_c), 'example', big_c[:3])
big_a = [ (b['E'], cg['ring']['aout'], cg['ring']['chand']) for b in blocks for cg in b['congs']
          if abs(cg['ring']['aout']) > 1 ]
print('congruences with |a|>1 outside:', len(big_a), 'example', big_a[:3])
per = Counter(sum(1 for cg in b['congs'] if abs(cg['ring']['chand'])>1) for b in blocks)
print('per-block count of |c|>1 congruences:', per.most_common())
