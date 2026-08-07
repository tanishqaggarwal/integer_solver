"""U30: the seven cheapest slots are exactly the ancestors of leaf 2081 (exp 235).
At each, with the cheapest lying leaf, enumerate EVERY honest leaf on the other side --
the slice through which the deliverable's 5-equation discount was found at the ROOT.
"""
import sys, time, collections, pickle
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentU_work')
import u20_sweep as S
import umodel as U

SP = pickle.load(open('u_slotprice.pkl', 'rb'))
cheap = sorted([b for b, s in SP.items() if s['min'] <= 19], key=lambda b: SP[b]['min'])
print('slots with generic minimum <= 19: %d  -> %s' % (len(cheap), cheap))
# is leaf 2081 under each?
for b in cheap:
    print('   beta=%-6d d=%-2d min=%-3d argmin_lying=%-6d  2081 in sub: %s'
          % (b, U.depth[b], SP[b]['min'], SP[b]['argmin_lying_leaf'], 2081 in U.sub[b]))

out = {}
t0 = time.time()
for b in cheap:
    lying = SP[b]['argmin_lying_leaf']
    la, lb = U.tree[b]
    side = lb if lying in U.LIVELEAF[la] else la      # honest leaves live on the OTHER side
    H = sorted(U.LIVELEAF[side])
    if b == U.ROOT:
        prof = pickle.load(open('u_honest_root.pkl', 'rb'))
    else:
        prof = {}
        for h in H:
            a, bb = (h, lying) if lying in U.LIVELEAF[lb] else (lying, h)
            prof[h] = S.price(b, a, bb, h)[0]
    out[b] = prof
    c = collections.Counter(prof.values())
    mn = min(prof.values())
    generic = c.most_common(1)[0][0]
    print('beta=%-6d d=%-2d lying=%-6d honest choices=%-4d  generic=%-3d  MIN=%-3d '
          'discount=%-2d  at honest %s   dist=%s  (%.0fs)'
          % (b, U.depth[b], lying, len(prof), generic, mn, generic - mn,
             [k for k, v in prof.items() if v == mn][:6], sorted(c.items())[:6], time.time() - t0))
    sys.stdout.flush()
pickle.dump(out, open('u_cheapslots.pkl', 'wb'))
allmin = min(min(p.values()) for p in out.values())
print('\nMINIMUM over every cheap-slot slice: %d' % allmin)
print('any evaluation below 7: %s' % any(v < 7 for p in out.values() for v in p.values()))
