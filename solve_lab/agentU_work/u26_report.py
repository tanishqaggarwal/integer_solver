"""U26: assemble the per-slot route price from the exhaustive shards."""
import sys, pickle, collections, json
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentU_work')
import umodel as U

rows = []
for i in range(4):
    rows += pickle.load(open('u_exact_%d.pkl' % i, 'rb'))
print('evaluations: %d  (expected %d)' % (len(rows), sum(
    len(U.LIVELEAF[U.tree[b][0]]) + len(U.LIVELEAF[U.tree[b][1]]) for b in U.SLOTS)))
bad = [r for r in rows if r[3] < 0]
print('errored: %d' % len(bad))

per = collections.defaultdict(list)
for beta, lying, honest, n in rows:
    if n >= 0:
        per[beta].append((n, lying))

slot = {}
for beta, v in per.items():
    v.sort()
    slot[beta] = {'min': v[0][0], 'argmin_lying_leaf': v[0][1], 'max': v[-1][0],
                  'nleaf': len(v), 'depth': U.depth[beta],
                  'nI': len(U.LIVELEAF[U.tree[beta][0]]),
                  'nJ': len(U.LIVELEAF[U.tree[beta][1]]),
                  'all': [x[0] for x in v]}
mins = [s['min'] for s in slot.values()]
print('\n=== PER-SLOT ROUTE PRICE (minimum over all leaf pairs and all carried values) ===')
print('slots priced: %d of %d merge slots' % (len(slot), len(U.SLOTS)))
print('min %d   median %d   mean %.1f   max %d' % (min(mins), sorted(mins)[len(mins)//2],
                                                   sum(mins)/len(mins), max(mins)))
print('distribution of per-slot minima:')
for k, c in sorted(collections.Counter(mins).items()):
    print('   %3d failing : %3d slots  %s' % (k, c, '#' * min(c, 60)))
print('\nslots at or below 10:')
for b, s in sorted(slot.items(), key=lambda kv: kv[1]['min'])[:20]:
    print('  beta=%-6d depth=%-2d |I|=%-3d |J|=%-3d  min=%-3d (lying leaf %d)  max=%d'
          % (b, s['depth'], s['nI'], s['nJ'], s['min'], s['argmin_lying_leaf'], s['max']))
print('\nROOT: %s' % slot.get(U.ROOT))
print('\nglobal minimum over EVERY evaluation: %d' % min(r[3] for r in rows if r[3] >= 0))
print('evaluations strictly below 7: %d' % sum(1 for r in rows if 0 <= r[3] < 7))
print('evaluations equal to 7      : %d' % sum(1 for r in rows if r[3] == 7))

# stratify by depth and by support size
print('\n=== stratified by slot depth ===')
bd = collections.defaultdict(list)
for b, s in slot.items():
    bd[s['depth']].append(s['min'])
for d in sorted(bd):
    v = sorted(bd[d])
    print('  depth %-2d : %3d slots  min %-3d median %-3d max %d' % (d, len(v), v[0], v[len(v)//2], v[-1]))
print('\n=== stratified by support size |I|+|J| ===')
bs = collections.defaultdict(list)
for b, s in slot.items():
    n = s['nI'] + s['nJ']
    k = 2 if n <= 2 else (4 if n <= 4 else (8 if n <= 8 else (16 if n <= 16 else (
        32 if n <= 32 else (64 if n <= 64 else 256)))))
    bs[k].append(s['min'])
for k in sorted(bs):
    v = sorted(bs[k])
    print('  |I|+|J| <= %-3d : %3d slots  min %-3d median %-3d max %d' % (k, len(v), v[0], v[len(v)//2], v[-1]))
pickle.dump(slot, open('u_slotprice.pkl', 'wb'))
json.dump({str(k): {kk: vv for kk, vv in s.items() if kk != 'all'} for k, s in slot.items()},
          open('u_slotprice.json', 'w'), indent=0)
print('\nwrote u_slotprice.pkl / u_slotprice.json')
