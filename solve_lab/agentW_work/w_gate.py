"""W stage 5: the liveness gate.  (1) is the gate the same var in all 3 congruences of a
block?  (2) is the SAME gate the one that multiplies the output pair (i5,i6) into the mux?
If yes, 'gate = 0' kills the law and the output propagation simultaneously -> not a family."""
import sys, os, re, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
import model
from collections import Counter
d = model.get(); A = d['atom_src']; AV = d['atom_vars']
byvar = {}
for i, vs in enumerate(AV):
    for v in vs: byvar.setdefault(v, []).append(i)
blocks = json.load(open('w_blocks3.json'))
def short(v): return [a for a in byvar.get(v, []) if len(A[a]) < 200]

same = Counter(); mux = Counter(); gatedef = Counter(); pin = Counter()
out_ok = 0; details = []
for b in blocks:
    Ls = set(cg['ring']['L'] for cg in b['congs'])
    same[len(Ls)] += 1
    L = next(iter(Ls)) if len(Ls) == 1 else None
    b['L'] = L
    if L is None: continue
    # how is L defined?
    ds = [A[a] for a in short(L) if re.fullmatch(r'x_%d - .*' % L, A[a])]
    kinds = set()
    for s in ds:
        if re.fullmatch(r'x_%d - x_\d+ \* x_\d+' % L, s): kinds.add('prod')
        elif re.fullmatch(r'x_%d - (0|1)' % L, s): kinds.add('pin' + s.split('- ')[1])
        elif re.fullmatch(r'x_%d - x_\d+' % L, s): kinds.add('alias')
        elif re.fullmatch(r'x_%d - \(1 - x_\d+\)' % L, s): kinds.add('not')
        else: kinds.add('other:' + s[:60])
    gatedef[tuple(sorted(kinds))] += 1
    # the output mux: consumers of i5 and i6 of the form  x_? - x_g * x_i5
    def gates_of(v):
        g = set()
        for a in short(v):
            m = re.fullmatch(r'x_(\d+) - x_(\d+) \* x_%d' % v, A[a])
            if m: g.add(int(m.group(2)))
            m = re.fullmatch(r'x_(\d+) - x_%d \* x_(\d+)' % v, A[a])
            if m: g.add(int(m.group(2)))
        return g
    g5, g6 = gates_of(b['i5']), gates_of(b['i6'])
    ok = (g5 == {L} and g6 == {L})
    out_ok += ok
    mux[(tuple(sorted(g5 - {L})), tuple(sorted(g6 - {L})), L in g5, L in g6)] += 1
    if not ok and len(details) < 5: details.append((b['E'], b['i5'], b['i6'], sorted(g5), sorted(g6), L))
print('blocks where all 3 congruences share one gate var:', dict(same))
print('gate definition kinds:', gatedef.most_common(8))
print('blocks where output (i5,i6) is multiplied by EXACTLY the gate L and nothing else:', out_ok, '/', len(blocks))
print('mux gate signature histogram:', mux.most_common(6))
for t in details: print('  exception:', t)
json.dump(blocks, open('w_blocks4.json', 'w'))
