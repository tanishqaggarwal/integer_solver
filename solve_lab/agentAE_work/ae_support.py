#!/usr/bin/env python3
"""agent AE -- the coordinator's item 3: selector support of the 927 c>1 lift conditions.

QUESTION.  Does U's partition theorem (no proper slot support has maskval >= N, so the fold
never wraps) make the 927 integer-lift conditions |S|-blind, killing agent AB's Sec.8
("an instance-side constraint on |S|") for all |S| at once?

METHOD.  A condition is a predicate on the wire values feeding it; a wire's value is a function
of the ON-set S only through S restricted to that wire's SELECTOR SUPPORT.  So measure, for each
of the 927, the selector support of the guard atom that carries the condition.

DEPENDENCIES, stated up front (I did not re-derive these from EQUATIONS.txt):
  * F's 39,033-atom parse            (agentT_work/mirror/F/circ4.pkl, circ2.vars_of)
  * L's handle/cofactor list          (agentT_work/mirror/L/handles.pkl)
  * U's per-wire selector-support closure (agentU_work/v_supp2.pkl)
Each is independently sanity-checked here before use; the join and every conclusion are mine.
"""
import os, sys, re, pickle, collections, json

LAB = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
FM = os.path.join(LAB, 'agentT_work', 'mirror', 'F')
LM = os.path.join(LAB, 'agentT_work', 'mirror', 'L')
sys.path.insert(0, FM)
from circ2 import vars_of                                    # noqa: E402

print('=' * 78)
print('STEP 1 -- U per-wire selector support: independent sanity checks')
U = pickle.load(open(os.path.join(LAB, 'agentU_work', 'v_supp2.pkl'), 'rb'))
supp = U['supp']
nonempty = {w: s for w, s in supp.items() if s}
fams = collections.Counter(frozenset(s) for s in nonempty.values())
distinct = sorted(fams, key=len)
allexp = set()
for s in distinct: allexp |= set(s)
print('  wires with a non-empty selector support : %d' % len(nonempty))
print('  distinct non-empty supports             : %d   (U claims 511 = 2*256-1)' % len(distinct))
print('  union of all supports                   : %d exponents, min %d max %d'
      % (len(allexp), min(allexp), max(allexp)))
big = [s for s in distinct if len(s) == 256]
print('  supports of size 256 (root)             : %d' % len(big))
sizes = sorted((len(s) for s in distinct), reverse=True)
print('  10 largest support sizes                : %s' % sizes[:10])
print('  singleton supports                      : %d   (expect 256, one per selector)'
      % sum(1 for s in distinct if len(s) == 1))
# laminarity: every pair nested or disjoint
viol = 0
ds = sorted(distinct, key=len)
for i in range(len(ds)):
    for j in range(i + 1, len(ds)):
        a, b = ds[i], ds[j]
        if a & b and not (a <= b):
            viol += 1
            if viol < 4: print('   LAMINARITY VIOLATION %s vs %s' % (sorted(a)[:6], sorted(b)[:6]))
print('  laminarity violations                   : %d   (U claims 0)' % viol)
if big:
    root = big[0]
    kids = [s for s in distinct if s < root and not any(s < t < root for t in distinct)]
    print('  maximal proper subsupports of the root   : %s' % sorted(len(k) for k in kids))

print()
print('=' * 78)
print('STEP 2 -- rebuild the 927 (agent T\'s identification, re-run here)')
d = pickle.load(open(os.path.join(FM, 'circ4.pkl'), 'rb'))
atoms = d['atoms']; names = list(atoms)
v2a = collections.defaultdict(list)
for i, a in enumerate(names):
    for u in vars_of(atoms[a]): v2a[u].append(i)
H = pickle.load(open(os.path.join(LM, 'handles.pkl'), 'rb'))
Us = sorted(set(H['handle']))
print('  cofactor wires u (from L\'s handle list) : %d' % len(Us))
defpat = re.compile(r'^\(x(\d+)-\(x(\d+)\*x(\d+)\)\)$')
rows = []; stats = collections.Counter()
for u in Us:
    ai = v2a[u]
    if len(ai) != 1: stats['u not in exactly 1 atom'] += 1; continue
    s = names[ai[0]].replace(' ', '')
    m = defpat.match(s)
    if not m: stats['definition not (h-(P*u))'] += 1; continue
    h = int(m.group(1)); Pv = int(m.group(2))
    others = [j for j in v2a[h] if j != ai[0]]
    if len(others) != 1: stats['h not in exactly 2 atoms'] += 1; continue
    g = names[others[0]].replace(' ', '')
    mm = re.search(r'\((\d+)\*x%d(?![0-9])\)' % h, g)
    if mm: c = int(mm.group(1))
    elif re.search(r'(?<![0-9])x%d(?![0-9])' % h, g): c = 1
    else: stats['h not found in guard'] += 1; continue
    stats['ok'] += 1
    rows.append(dict(u=u, h=h, P=Pv, c=c, gi=others[0], di=ai[0]))
c1 = sum(1 for r in rows if r['c'] == 1); cg = [r for r in rows if r['c'] > 1]
print('  parse outcome: %s' % dict(stats))
print('  c == 1 : %d      c > 1 : %d   (L/P/T report 2754 / 927)' % (c1, len(cg)))

print()
print('=' * 78)
print('STEP 3 -- selector support of each c>1 condition')
def wires_support(ws):
    s = set(); miss = 0
    for w in ws:
        if w in supp: s |= supp[w]
        else: miss += 1
    return s, miss

def guard_support(r):
    ws = vars_of(atoms[names[r['gi']]]) | vars_of(atoms[names[r['di']]])
    return wires_support(ws)

hist = collections.Counter(); missing = 0; per = []
for r in cg:
    s, miss = guard_support(r)
    missing += miss
    hist[len(s)] += 1
    per.append((r['u'], r['h'], r['c'], len(s), frozenset(s)))
print('  wires in a guard with no entry in U\'s supp table: %d' % missing)
print('  support-size histogram over the %d c>1 conditions:' % len(cg))
for k in sorted(hist): print('     |support| = %3d : %5d conditions' % (k, hist[k]))
mx = max(h for h in hist)
print('  MAX selector support over the 927 : %d' % mx)
print('  conditions with support = all 256 : %d' % hist.get(256, 0))
distinct_supports = collections.Counter(x[4] for x in per)
print('  distinct supports carrying a c>1 condition : %d' % len(distinct_supports))

json.dump(dict(hist={str(k): v for k, v in hist.items()},
               n_cgt1=len(cg), n_c1=c1, max_support=mx,
               n_root_support=hist.get(256, 0),
               n_distinct_supports=len(distinct_supports),
               u_distinct_supports=len(distinct), u_laminarity_violations=viol),
          open('res_support.json', 'w'), indent=1)
print('\nwrote res_support.json')
