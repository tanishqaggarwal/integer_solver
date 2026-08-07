import sys, collections, os
from model import Model, load_assign
HERE = os.path.dirname(os.path.abspath(__file__))
M = Model()
v = load_assign(os.path.join(HERE, '..', 'best', 'new_instance_partial_39026.json'))
P = 2**256 - 2**32 - 977
var2atoms = collections.defaultdict(list)
for i, vs in enumerate(M.avars):
    for x in vs:
        var2atoms[x].append(i)

roots = [int(a) for a in sys.argv[1:]]
seen = set()
frontier = list(roots)
depth = 0
MAXD = int(os.environ.get('MAXD', '3'))
while frontier and depth <= MAXD:
    nxt = []
    print(f"===== depth {depth}: {len(frontier)} vars")
    for x in frontier:
        if x in seen:
            continue
        seen.add(x)
        print(f"  X{x} = {v[x]}   (mod p = {v[x] % P})")
        for a in var2atoms[x]:
            print(f"      a{a}: {M.src[a][:170]}")
            for y in M.avars[a]:
                if y not in seen:
                    nxt.append(y)
    frontier = nxt
    depth += 1
