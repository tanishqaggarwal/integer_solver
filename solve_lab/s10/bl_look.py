import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256 - 2**32 - 977
definer, atom_out = L.definer, L.atom_out
FREE=set(ad.FREE)

def show(a, maxlen=400):
    src = L.atom_src[a]
    out = atom_out.get(a)
    print(f'a{a}: out={out} nmon={len(L.polys[a])} vars={sorted(L.avars[a])[:14]}')
    print(f'    src: {src[:maxlen]}')

for a in [21617, 29539, 22229, 22230, 35758, 35759, 35760, 35761, 35762, 7930, 22231,
          37662, 40826]:
    show(a)
    print()

print('=== MUX region ===')
for t in (21279, 7075, 2081, 4287, 7068, 2099):
    d = definer.get(t)
    print(f'x_{t}: definer a{d}' + (f'  src={L.atom_src[d][:220]}' if d is not None else '  FREE'))

print('\n=== shape census of all atoms ===')
shape = collections.Counter()
for a, p in enumerate(L.polys):
    deg = max((len(m) for m in p), key=lambda x: x) if p else 0
    shape[(len(p), max((len(m) for m in p), default=0))] += 1
for k, c in sorted(shape.items(), key=lambda kv: -kv[1])[:25]:
    print(f'  (nmon={k[0]}, maxdeg={k[1]}): {c}')

print('\n=== atoms containing a big constant (|c|>2^32) as a bare monomial ===')
bigc = []
for a, p in enumerate(L.polys):
    for m, c in p.items():
        if len(m) == 0 and abs(c) > 2**32:
            bigc.append((a, c)); break
print(f'  {len(bigc)} atoms')
for a, c in bigc[:12]:
    print(f'   a{a}: nmon={len(L.polys[a])} const={str(c)[:40]}... src={L.atom_src[a][:150]}')
