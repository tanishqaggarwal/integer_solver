"""S11 step 1: what IS this instance? Bit-length census, constant identification,
atom taxonomy."""
import os, sys, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, ad
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)
v = L.load(os.path.join(HERE, 'mod9118_0.json'))

print('=== bit-length census ===')
def hist(xs, label):
    h = collections.Counter()
    for x in xs:
        b = x.bit_length()
        h[0 if b == 0 else (1 if b == 1 else (256 if b <= 256 else
            (300 if b <= 300 else (600 if b <= 600 else
             (1200 if b <= 1200 else 3000)))))] += 1
    print(f'  {label}: {dict(sorted(h.items()))}')
hist([v[u] for u in range(L.NVARS)], 'all variables')
hist([v[u] for u in FREE], 'free inputs   ')
hist([v[u] for u in range(L.NVARS) if u not in FREE], 'gate outputs  ')
big = sorted((u for u in FREE if v[u].bit_length() > 600),
             key=lambda u: -v[u].bit_length())
print(f'  free inputs over 600 bits: {len(big)} -> {big[:20]}')

print('\n=== atom taxonomy ===')
tax = collections.Counter()
for a in range(L.NA):
    poly = L.polys[a]
    ne = len(L.atom2eq[a])
    ms = list(poly)
    if len(ms) == 2 and any(len(m) == 2 and m[0] == m[1] for m in ms): tax['boolean'] += 1
    elif len(ms) == 3 and all(len(m) == 1 for m in ms): tax['linear gadget'] += 1
    elif ne == 1: tax['1-equation check'] += 1
    elif len(set(L.avars[a])) <= 3: tax['small (<=3 vars)'] += 1
    else: tax['other'] += 1
print(f'  {dict(tax)}')
print(f'  gate atoms {len(atom_out)}, check atoms {L.NA - len(atom_out)}')
ne = collections.Counter(len(L.atom2eq[a]) for a in range(L.NA))
print(f'  atoms by #equations: {dict(sorted(ne.items()))}')

print('\n=== broadcast constants ===')
cnt = collections.Counter(v)
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
for val, n in cnt.most_common(6):
    if n < 40 or val in (0, 1): continue
    print(f'  x{n}: {val}')
    print(f'      hex {hex(val)}')
    print(f'      bits {val.bit_length()}  mod p {val % P}')
    print(f'      mod n {val % N}   == n? {val == N}')
    if val.bit_length() > 256:
        q, r = divmod(val, P)
        print(f'      = {q}*p + {r}')
        print(f'      q hex {hex(q)}  r hex {hex(r)}')
