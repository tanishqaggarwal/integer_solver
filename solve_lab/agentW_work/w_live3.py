"""W stage 6c: the liveness cone is a BOOLEAN circuit.  Recognise AND/OR/NOT/alias/pin
structurally (OR = (a+b) - a*b on the SAME pair) and prove every gate lands in {0,1}."""
import sys, os, re, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
import model
from collections import Counter, deque
d = model.get(); A = d['atom_src']; AV = d['atom_vars']
byvar = {}
for i, vs in enumerate(AV):
    for v in vs: byvar.setdefault(v, []).append(i)
blocks = json.load(open('w_blocks4.json'))
def short(v): return [a for a in byvar.get(v, []) if len(A[a]) < 200]
def defs(v): return [A[a] for a in short(v) if re.fullmatch(r'x_%d - .*' % v, A[a])]
BOOL = set()
for s in A:
    m = re.fullmatch(r'x_(\d+) \* x_(\d+) - x_(\d+)', s)
    if m and m.group(1) == m.group(2) == m.group(3): BOOL.add(int(m.group(1)))
    m = re.fullmatch(r'x_(\d+) \* \(x_(\d+) - 1\)', s)
    if m and m.group(1) == m.group(2): BOOL.add(int(m.group(1)))
    m = re.fullmatch(r'2 \* x_(\d+) \* \(1 - x_(\d+)\)', s)
    if m and m.group(1) == m.group(2): BOOL.add(int(m.group(1)))
RULES = [('prod', r'x_%d - x_(\d+) \* x_(\d+)'), ('not', r'x_%d - \(1 - x_(\d+)\)'),
         ('alias', r'x_%d - x_(\d+)'), ('sum', r'x_%d - \(x_(\d+) \+ x_(\d+)\)'),
         ('diff', r'x_%d - \(x_(\d+) - x_(\d+)\)'), ('pin', r'x_%d - (-?\d+)')]
node = {}
gates = set(b['L'] for b in blocks)
Q = deque(gates); seen = set()
while Q:
    v = Q.popleft()
    if v in seen: continue
    seen.add(v)
    got = None
    for k, pat in RULES:
        for s in defs(v):
            m = re.fullmatch(pat % v, s)
            if m:
                got = (k, [int(g) for g in m.groups()]); break
        if got: break
    node[v] = got or ('leaf', [])
    if got and got[0] != 'pin': Q += got[1]
# classify every node as a boolean gate
kind = {}; unresolved = []
for v in seen:
    k, ar = node[v]
    if k == 'leaf': kind[v] = ('LEAF', ar) if v in BOOL else ('FREE', ar)
    elif k == 'pin': kind[v] = ('CONST', ar)
    elif k == 'alias': kind[v] = ('ALIAS', ar)
    elif k == 'prod': kind[v] = ('AND', ar)
    elif k == 'not': kind[v] = ('NOT', ar)
    elif k == 'diff':
        p, q = ar                                   # v = p - q
        kp, ap = node.get(p, ('?', [])); kq, aq = node.get(q, ('?', []))
        if kp == 'sum' and kq == 'prod' and sorted(ap) == sorted(aq): kind[v] = ('OR', sorted(ap))
        else: unresolved.append((v, 'diff', A[short(v)[0]])); kind[v] = ('?', ar)
    elif k == 'sum': kind[v] = ('SUM', ar)
    else: kind[v] = ('?', ar)
print('cone size', len(seen), Counter(kind[v][0] for v in seen).most_common())
print('unresolved diff nodes:', len(unresolved), unresolved[:3])
# are SUM / raw-prod nodes used ONLY inside an OR pattern?
orparts = set()
for v in seen:
    if kind[v][0] == 'OR':
        p, q = node[v][1]; orparts.add(p); orparts.add(q)
sums = [v for v in seen if kind[v][0] == 'SUM']
print('SUM nodes:', len(sums), ' all inside an OR:', all(v in orparts for v in sums))
# consumers of SUM nodes outside the cone?
esc = [v for v in sums if any(a for a in byvar[v] if len(A[a]) < 200 and not re.fullmatch(r'x_%d - .*'%v, A[a]) and v not in orparts)]
# boolean closure
BOOLOK = set(v for v in seen if kind[v][0] == 'LEAF')
BOOLOK |= set(v for v in seen if kind[v][0] == 'CONST' and kind[v][1][0] in (0, 1))
for _ in range(200):
    add = set()
    for v in seen:
        if v in BOOLOK: continue
        k, ar = kind[v]
        if k in ('ALIAS', 'NOT') and ar[0] in BOOLOK: add.add(v)
        elif k in ('AND', 'OR') and all(p in BOOLOK for p in ar): add.add(v)
    if not add: break
    BOOLOK |= add
print('nodes proved in {0,1}:', len(BOOLOK), 'of', len(seen))
notb = sorted(set(seen) - BOOLOK)
print('not proved boolean:', len(notb), Counter(kind[v][0] for v in notb).most_common())
gb = [b['E'] for b in blocks if b['L'] not in BOOLOK]
print('GATES not proved boolean:', len(gb), gb[:5])
print('gate kinds:', Counter(kind[b['L']][0] for b in blocks).most_common())
