"""Trade-knob walk: generate valid test cases BY CLASS instead of sampling classes and hoping.

Base: cfg0 near-solution (cfg0 + the particular solution n0; bad = {a20215, a28647}, 39,005).
Moves: the 1-for-1 trade knobs x_14853 (a20212 -1 / a28647 +1), x_6083 (a7389 +1 / a28647 -1),
x_31339 (a10187 +1 / a20215 -1), x_18956 (a747 +8863713 / a20215 +1).  After each move the OTHER
rows are re-solved exactly by lat3.analyse, which then tests membership on the targets.

REPORTED PRIMARY: the count of VALID cases (other rows re-solvable AND the mod-p class actually
moved).  The starvation rate is what decides whether the question is answerable this way.

FIRST it checks the thing that would make this whole experiment vacuous: whether the trade knobs
already lie inside the 54-knob affine span that lat3.analyse optimises over.  If they do, moving
along them cannot change the membership answer and the walk proves nothing -- the same defect
Result B had, and better to find it up front than to report it afterwards.
"""
import sys, json, collections, pickle, time, random
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentS_work')
import common as C, lat2, lat3
import harness as H, engine as E, fast, intsolve
P = C.P
TRADE = [14853, 6083, 31339, 18956]
TGT = [20215, 28647]

print("=== step 0: are the trade knobs inside the affine span lat3 already optimises over? ===", flush=True)
seed0 = dict(C.BASE)
v0, bad0, aff, atoms, hs = lat2.system(seed0)
knobs = sorted(aff)
inside = [f for f in TRADE if f in aff]
print("54-knob affine set size = %d ; trade knobs present in it: %s ; absent: %s"
      % (len(knobs), inside, [f for f in TRADE if f not in aff]), flush=True)
if len(inside) == len(TRADE):
    print("  ALL FOUR trade knobs are already inside the span.", flush=True)
    print("  => lat3.analyse already explores every integer combination of them, so displacing", flush=True)
    print("     along a trade knob and re-solving CANNOT change the membership answer.", flush=True)
    print("  Running the walk anyway to confirm empirically rather than argue it.", flush=True)

res = lat3.analyse(seed0, 'cfg0')
y, ker, n0, knobs, aff, v0, bad0, atoms = res
base = dict(seed0)
for j, f in enumerate(knobs):
    if n0[j]:
        base[f] = v0[f] + n0[j]
vb = E.forward(base)
badb = E.badatoms(vb)
nsb = {'v': vb, '__builtins__': {}}
class0 = tuple(eval(H.acodes[a], nsb) % P for a in TGT)
print("\nbase near-solution: bad=%s ; class=(%s..., %s...)"
      % (sorted(badb), str(class0[0])[:22], str(class0[1])[:22]), flush=True)

random.seed(77)
moves = []
for f in TRADE:
    for n in [1, 2, 3, -1, -5, 12345, random.randint(1, 10**9), random.randint(1, 10**40)]:
        moves.append(({f: n}, 'x_%d += %s' % (f, str(n)[:14])))
for t in range(8):
    m = {f: random.randint(-10**6, 10**6) for f in random.sample(TRADE, 2)}
    moves.append((m, 'pair ' + '+'.join('x_%d' % f for f in m)))

valid = blocked = solved = 0
notnear = 0
sameclass = 0
classes = set()
print("\n=== walk: %d moves ===" % len(moves), flush=True)
for mv, tag in moves:
    ns = dict(base)
    for f, n in mv.items():
        ns[f] = ns.get(f, v0[f]) + n
    try:
        v = E.forward(ns); bad = E.badatoms(v)
        nsx = {'v': v, '__builtins__': {}}
        cls = tuple(eval(H.acodes[a], nsx) % P for a in TGT)
    except Exception as e:
        print("  [%s] forward ERR" % tag, flush=True); continue
    moved = (cls != class0)
    classes.add(cls)
    try:
        r = lat3.analyse(ns, '  ' + tag)
    except Exception as e:
        print("  [%s] analyse ERR %s" % (tag, type(e).__name__), flush=True); continue
    if r is None:
        notnear += 1
        print("  [%s] class-moved=%s -> OTHER ROWS INFEASIBLE (not a valid case)" % (tag, moved), flush=True)
        continue
    if not moved:
        sameclass += 1
        print("  [%s] other rows solvable but CLASS UNCHANGED (not a valid case)" % tag, flush=True)
        continue
    valid += 1
    if r[0] is not None:
        solved += 1
        print("  [%s] *** VALID and SOLVED -- obstruction is cfg0-local ***" % tag, flush=True)
        json.dump({str(x): str(int(z)) for x, z in ns.items()}, open('S_trade_hit.json', 'w'))
    else:
        blocked += 1
        print("  [%s] VALID, blocked" % tag, flush=True)

print("\n=== SUMMARY (trade walk) ===", flush=True)
print("moves attempted            : %d" % len(moves), flush=True)
print("other rows infeasible      : %d" % notnear, flush=True)
print("solvable but class unmoved : %d" % sameclass, flush=True)
print("VALID CASES                : %d   (blocked=%d solved=%d)" % (valid, blocked, solved), flush=True)
print("distinct target classes hit: %d" % len(classes), flush=True)
if valid == 0:
    print("=> the trade knobs STARVE TOO. This line is closed: it cannot answer the question.", flush=True)
elif solved == 0:
    print("=> %d valid cases, all blocked." % valid, flush=True)
