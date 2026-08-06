"""Extract the full LOAD TABLE: for every boolean free input, the constants it gates in."""
import sys, os, json, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
P = L.P
FREE = set(u for u in range(L.NVARS) if L.definer.get(u) is None)

# a load pin looks like:  bit*X  -  HUGE*bit  -  c*handle      (bit, X both free)
loads = collections.defaultdict(list)   # bit -> [(atom, X, HUGE, coeff_of_handle, handle)]
for a in range(L.NA):
    if L.atom_out.get(a) is not None:
        continue
    Pp = L.polys[a]
    # find a degree-2 monomial of two free vars
    quads = [(mm, c) for mm, c in Pp.items() if len(mm) == 2 and mm[0] in FREE and mm[1] in FREE]
    if len(quads) != 1:
        continue
    mm, cq = quads[0]
    lin = {m[0]: c for m, c in Pp.items() if len(m) == 1}
    # the bit is the one that also appears linearly with a HUGE coefficient
    for bit, other in ((mm[0], mm[1]), (mm[1], mm[0])):
        if bit in lin and abs(lin[bit]) > 10**30:
            loads[bit].append((a, other, -lin[bit] // cq if cq else None, lin[bit], cq))
            break

print(f"boolean free inputs carrying load pins: {len(loads)}")
tot = sum(len(v) for v in loads.values())
print(f"total load pins: {tot}")
cnt = collections.Counter(len(v) for v in loads.values())
print("pins per bit:", dict(cnt))

# collect the constants
consts = {}
for b, lst in loads.items():
    consts[b] = sorted(set(abs(h) for _, _, _, h, _ in lst))
allc = sorted(set(c for v in consts.values() for c in v))
print(f"distinct load constants: {len(allc)}")
print("bit-count of constants:", collections.Counter(c.bit_length() for c in allc[:50]))
res = [c % P for c in allc]
print(f"constants >= P: {sum(1 for c in allc if c >= P)} of {len(allc)}")
print(f"distinct residues mod P: {len(set(res))}")
json.dump({'loads': {str(k): [[a, o, None, h, cq] for a, o, _, h, cq in v] for k, v in loads.items()}},
          open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'loads.json'), 'w'))
print("saved data/loads.json")
