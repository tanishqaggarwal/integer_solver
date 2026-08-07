"""U32: anatomy of the 5-equation discount {2554,6816,8124,8680,9421} and of the 11 DRV knobs."""
import sys, collections, pickle
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentU_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentU_work/mirror')
import harness as H
import engine3 as E3
import umodel as U

ENG = E3.Eng(E3.BASE_DEMOTE)
FREE = set(ENG.FREE)
DISC = [2554, 6816, 8124, 8680, 9421]
DRV = [642, 1329, 8731, 9118, 9413, 10903, 17325, 18956, 28730, 29854, 31864]

print('=== the five discount equations ===')
touch = collections.defaultdict(set)
for e in DISC:
    issq, outer, terms = H.eqt[e]
    ats = [a for c, a in terms if a >= 0]
    vs = set()
    for a in ats:
        vs |= set(H.avars[a])
    fv = sorted(vs & FREE)
    print(' eq %-5d sq=%s atoms=%-3d vars=%-4d free=%-4d  DRV in it: %s'
          % (e, issq, len(ats), len(vs), len(fv), sorted(set(DRV) & vs)))
    for d in set(DRV) & vs:
        touch[d].add(e)
print('\nDRV variable -> which discount equations it appears in:')
for d in DRV:
    print('  x_%-6d free=%-5s occurs in %-2d atoms ; discount eqs %s'
          % (d, d in FREE, len(H.occ.get(d, [])), sorted(touch.get(d, []))))

print('\n=== atoms of the discount equations that contain a DRV variable ===')
seen = set()
for e in DISC:
    issq, outer, terms = H.eqt[e]
    for c, a in terms:
        if a < 0 or a in seen:
            continue
        vs = set(H.avars[a])
        if vs & set(DRV):
            seen.add(a)
            s = H.atoms[a]
            print('  eq%-5d atom %-6d free-in-atom=%s' % (e, a, sorted(vs & FREE)))
            print('        %s' % (s if len(s) < 260 else s[:260] + '...'))

# which atoms can each DRV var be solved from?
print('\n=== solvable-from atoms per DRV variable (atom, is-it-a-definer-elsewhere) ===')
for d in DRV:
    rows = []
    for i in H.occ.get(d, []):
        s = H.atoms[i]
        rows.append((i, len(s), E3.ATOM2VAR.get(i)))
    print('  x_%-6d %s' % (d, rows))
pickle.dump({'disc': DISC, 'drv': DRV, 'touch': dict(touch)}, open('u_drvmap.pkl', 'wb'))
