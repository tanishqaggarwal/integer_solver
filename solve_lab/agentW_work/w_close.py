"""TASK 1d: close the frame-B budget at every j, using the enumerated cocircuits.

Lemma (proved, not assumed).  Let x be any knob delta and T = N(x) the set of SAT rows it
breaks.  SAT rows are homogeneous, so x lies in ker_Z(A_SAT\T) = Z^34 cap ker_Q(A_SAT\T).
Let T' be the union of the minimal cocircuits contained in T.  If rank(A_SAT\T) equalled
rank(A_SAT\T') the two kernels would coincide and x would satisfy the rows of T\T', i.e.
T = T'.  Otherwise T\T' contains a cocircuit C' of the deletion minor M\T', and every such
C' is C cap (E\T') for a cocircuit C of M with C contained in T -- hence C contained in T',
hence C' empty.  Contradiction.  So **T is always a union of minimal cocircuits**, and it is
enough to test unions of minimal cocircuits of total size <= 6.
"""
import sys, os, itertools, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w_setup as S

cc = [frozenset(int(x) if x != 'S' else 'S' for x in T) for T in json.load(open('w_cocirc_raw.json'))]
minimal = [c for c in cc if not any(d < c for d in cc)]
minimal = sorted(set(minimal), key=lambda c: (len(c), sorted(map(str, c))))
print('cocircuit-containing sets enumerated: %d ; MINIMAL ones: %d' % (len(cc), len(minimal)))
from collections import Counter as _C
print('  minimal by size:', dict(sorted(_C(len(c) for c in minimal).items())))
for c in minimal:
    if len(c) <= 2: print('   size %d : %s' % (len(c), sorted(map(str, c))))

# all unions of minimal cocircuits with total size <= 6
Bs = {frozenset()}
frontier = {frozenset()}
while frontier:
    nxt = set()
    for b in frontier:
        for c in minimal:
            u = b | c
            if len(u) <= 6 and u not in Bs:
                Bs.add(u); nxt.add(u)
    frontier = nxt
Bs = sorted(Bs, key=lambda b: (len(b), sorted(map(str, b))))
CAP = int(os.environ.get('WCLOSECAP', '400000'))
if len(Bs) > CAP:
    print('!! union closure capped at %d of %d' % (CAP, len(Bs))); Bs = Bs[:CAP]
print('\ncandidate break-sets (unions of minimal cocircuits, size <= 6): %d' % len(Bs))
from collections import Counter
print('  by size:', dict(sorted(Counter(len(b) for b in Bs).items())))

t0 = time.time(); res = {}; best = (0, None, None); n = 0
for k in range(1, 8):
    for P in itertools.combinations(S.FAIL, k):
        mb, mbB = None, None
        for b in Bs:
            if len(b) >= k: break            # cannot gain
            keep = [e for e in S.SAT if e not in b]
            n += 1
            if S.solve(keep + list(P)) is not None:
                mb, mbB = len(b), sorted(b, key=str); break
        gain = (k - mb) if mb is not None else None
        res[str(list(P))] = {'minbreak_lt_k': mb, 'B': [str(x) for x in mbB] if mbB else None, 'gain': gain}
        if gain and gain > best[0]:
            best = (gain, list(P), mbB)
            sol = S.solve([e for e in S.SAT if e not in mbB] + list(P))
            S.price(sol, 'buy %s break %s' % (list(P), mbB), tagfile='close')
        print('  buy %-40s : any break-set of size < %d feasible? %s'
              % (list(P), k, 'YES gain %d with B=%s' % (gain, [str(x) for x in mbB]) if gain else 'NO'), flush=True)
print('\n%d integer solves in %.0fs' % (n, time.time() - t0))
print('BEST GAIN over every bought-set and every rank-dropping break-set of size <= 6: %s' % (best,))
json.dump({'nminimal': len(minimal), 'nBs': len(Bs), 'res': res}, open('w_close.json', 'w'), indent=1)
