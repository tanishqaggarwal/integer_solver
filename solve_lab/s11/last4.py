"""At the closehit2 state: residuals a14445, a27139, a34580, a33796 oscillate.
   Scan every free input for controls, excluding the locked (gamma/mirror/arithmetic) set."""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
M = 8640431
LOCK = {31339, 33708, 490, 91, 19750, 7497, 22820, 14853, 14393, 11436, 14515, 16742,
        22162, 30213, 8386, 21868, 16441, 28955, 2751, 18751}
NAMES = ['a14445', 'a27139', 'a34580', 'a33796', 'mirror3719', 'mirror25118', 'gamma']

v0 = [0] * L.NVARS
for k, x in json.load(open(os.path.join(HERE, 'data', 'closehit2.json'))).items():
    v0[int(k)] = int(x)
fw.forward(v0)
print("bad:", fw.bad_checks(v0), "failing:", len(L.failing_eqs(L.all_atom_values(v0))))


def res(v):
    g = (v[12000] // P) % M if v[12000] % P == 0 else -1
    return [(v[33129] - v[3757]) % P, (v[37088] - v[13585]) % P,
            (v[33708] - v[10170]) % P, (v[31339] - v[6858]) % P,
            v[3719] % P, v[25118] % P, g]


base = res(v0)
print("residuals:", [('0' if x == 0 else 'nz') for x in base])
FREE = [u for u in range(L.NVARS) if L.definer.get(u) is None and u not in LOCK]
t0 = time.time()
rows = {i: [] for i in range(7)}
v = v0[:]
for n, u in enumerate(FREE):
    old = v[u]
    v[u] = old + 1
    fw.forward(v)
    r1 = res(v)
    v[u] = old
    fw.forward(v)
    for i in range(7):
        if r1[i] != base[i]:
            rows[i].append(u)
    if n % 1500 == 0:
        print(f"  {n}/{len(FREE)} ({time.time()-t0:.0f}s)", flush=True)
print(f"scan done ({time.time()-t0:.0f}s)")
for i in range(7):
    print(f"  {NAMES[i]:12s}: {len(rows[i])} -> {rows[i][:14]}")
json.dump({NAMES[i]: rows[i] for i in range(7)}, open(os.path.join(HERE, 'data', 'last4.json'), 'w'))
