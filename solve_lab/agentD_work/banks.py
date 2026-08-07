"""Split the 256 gating bits into the two OR-tree banks, then search for a pair
(i in bank1, j in bank2) satisfying the two addition identities exactly mod p."""
import json, sys, itertools
import dlib as L
import ortree
import io, contextlib
P = L.P

v = L.load('D_state1.json')
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    l1 = ortree.expand(7715, v)
    l2 = ortree.expand(34554, v)
t = json.load(open('table.json'))
bits = set(int(b) for b in t)
b1 = [u for u in l1 if u in bits]
b2 = [u for u in l2 if u in bits]
print('tree x_7715 leaves', len(l1), 'gating bits', len(b1))
print('tree x_34554 leaves', len(l2), 'gating bits', len(b2))
print('24601 in bank1', 24601 in b1, ' in bank2', 24601 in b2)
print('2081  in bank1', 2081 in b1, ' in bank2', 2081 in b2)

CB = 125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626 % P
CC = 91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002 % P
CA = 97553848499418123410591666447050222001188385549510401465815187079080512838891

def consts(b):
    return [C % P for a, x, C in t[str(b)]]

hits = []
for i in b1:
    ci = consts(i)
    for oi in (0, 1):
        X1, Y1 = ci[oi], ci[1 - oi]
        for j in b2:
            cj = consts(j)
            for oj in (0, 1):
                X2, Y2 = cj[oj], cj[1 - oj]
                A = ((X1 + X2 + CC + CA) * (X2 - X1) ** 2 - (Y2 - Y1) ** 2) % P
                A0 = ((X1 + X2 + CC) * (X2 - X1) ** 2 - (Y2 - Y1) ** 2) % P
                B = ((CB + Y1) * (X2 - X1) - (Y2 - Y1) * (X1 - CC)) % P
                if A == 0 or A0 == 0 or B == 0:
                    hits.append((i, oi, j, oj, A == 0, A0 == 0, B == 0))
print('hits:', len(hits))
for h in hits[:40]:
    print('  ', h)
json.dump({'bank1': b1, 'bank2': b2}, open('banks.json', 'w'))
