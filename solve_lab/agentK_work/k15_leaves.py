#!/usr/bin/env python3
"""K15: find every atom carrying a large integer literal -> the leaf constants and the target."""
import sys, os, json, re, collections
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
from cascadep import CascadeP, NV, P
from parse import node_str

C = CascadeP()
big = collections.defaultdict(list)
shapes = collections.Counter()
for i, nm in enumerate(C.names):
    for m in re.findall(r'(?<![x\d])(\d{20,})', nm):
        big[int(m)].append(i)
    if re.search(r'(?<![x\d])\d{20,}', nm):
        shapes[re.sub(r'(?<![x\d])\d{20,}', 'BIG', re.sub(r'x\d+', 'X', re.sub(r'(?<![x\d])\d{1,19}(?![\d])', 'c', nm)))] += 1
print('distinct large literals:', len(big))
print('shape histogram:')
for s, c in shapes.most_common(20): print('   %4d  %s' % (c, s))
vals = sorted(big)
print('literal magnitudes: min %d digits, max %d digits' % (len(str(vals[0])), len(str(vals[-1]))))
cnt = collections.Counter(len(str(x)) for x in vals)
print('digit-length histogram:', sorted(cnt.items()))
json.dump({str(k): v for k, v in big.items()}, open(K + '/bigliterals.json', 'w'))
# print a few example atoms per shape
byshape = collections.defaultdict(list)
for i, nm in enumerate(C.names):
    if re.search(r'(?<![x\d])\d{20,}', nm):
        s = re.sub(r'(?<![x\d])\d{20,}', 'BIG', re.sub(r'x\d+', 'X', re.sub(r'(?<![x\d])\d{1,19}(?![\d])', 'c', nm)))
        byshape[s].append(nm)
for s in list(byshape)[:12]:
    print('---', s, 'n=', len(byshape[s]))
    print('    ', byshape[s][0][:160])
