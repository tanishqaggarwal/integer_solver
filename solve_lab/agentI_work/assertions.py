"""Enumerate the mod-p assertion gadgets and the conditional-pin gadgets."""
import pickle, os, collections
from model import Model
HERE = os.path.dirname(os.path.abspath(__file__))
M = Model()
val = pickle.load(open(os.path.join(HERE, 'prop0.pkl'), 'rb'))
P = 2**256 - 2**32 - 977
pvars = set(v for v in range(len(val)) if val[v] == P)
zeros = set(v for v in range(len(val)) if val[v] == 0)
ones = set(v for v in range(len(val)) if val[v] == 1)
print("p-vars:", len(pvars), " zeros:", len(zeros), " ones:", len(ones))

polys = M.polys
assertions = []   # atoms with a monomial (pvar, h)
condpins = []
other = []
for i, q in enumerate(polys):
    has_p = False
    for m in q:
        if any(x in pvars for x in m):
            has_p = True
    if has_p:
        assertions.append(i)
print("atoms involving a p-valued variable:", len(assertions))
import re
sh = collections.Counter()
for i in assertions:
    t = re.sub(r'X\d+', 'V', M.src[i]); t = re.sub(r'\d+', 'N', t)
    sh[t] += 1
for k, n in sh.most_common(20):
    print(f"   {n:5d}  {k}")

# conditional pins: atoms with a huge literal constant
huge = []
for i, q in enumerate(polys):
    if any(abs(c) > 10**60 for c in q.values()):
        huge.append(i)
print("\natoms with a >2^200 literal:", len(huge))
sh = collections.Counter()
for i in huge:
    t = re.sub(r'X\d+', 'V', M.src[i]); t = re.sub(r'\d+', 'N', t)
    sh[t] += 1
for k, n in sh.most_common(20):
    print(f"   {n:5d}  {k}")
pickle.dump({'pvars': pvars, 'assertions': assertions, 'huge': huge},
            open(os.path.join(HERE, 'assert.pkl'), 'wb'))
