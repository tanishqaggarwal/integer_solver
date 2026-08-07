"""A negative result from a linear model is only sound if the model is exact.

My columns were validated with a two-point probe (+1, +2), which certifies affineness only
against a quadratic alternative.  Re-probe every (knob, reachable check) pair at t = 0,1,2,3,5,7
and require the value to be exactly affine in t.  Anything that fails here would have been
silently mismodelled.
"""
import sys, os, json, re
HERE = '/home/user/integer_solver/solve_lab/agentH_work'
sys.path.insert(0, HERE)
os.chdir(HERE)
sys.set_int_max_str_digits(20_000_000)
OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/fb_lincheck.log', 'w', buffering=1)


def say(*a):
    print(*a, file=LOG)


import frameB
VAR_RE = re.compile(r'x_(\d+)')
fr = frameB.Frame([642, 28730, 29854, 31864])
W = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
vw = [0] * frameB.NV
for k, val in W.items():
    vw[int(k[2:]) if k.startswith('x_') else int(k)] = int(val)
fv0 = {u: vw[u] for u in fr.free if vw[u] != 0}
st0 = frameB.State(fr, fv0)
A = 37887
src = frameB.atom_src[A]
half = (len(src) - 3) // 2
SCODE = compile(VAR_RE.sub(r'v[\1]', src[:half]), '<S>', 'eval')
Sval = lambda st: eval(SCODE, {'v': st.v, '__builtins__': {}})

NZ = sorted(st0.nz())
U = set()
for a in NZ:
    U |= set(fr.SUPV.get(a, []))
KNOB = sorted(U | set(fr.SUPV.get(A, [])))
base = {a: st0.av[a] for a in fr.checks}
S0 = Sval(st0)
TS = [1, 2, 3, 5, 7]

bad_nonaffine = {}
Sbad = []
for u in KNOB:
    o = fv0.get(u, 0)
    sts = {t: st0.clone().set_free({u: o + t}) for t in TS}
    d1 = {a: sts[1].av[a] - base[a] for a in fr.checks}
    for a in fr.checks:
        for t in TS:
            if sts[t].av[a] - base[a] != t * d1[a]:
                bad_nonaffine.setdefault(u, set()).add(a)
                break
    e1 = Sval(sts[1]) - S0
    for t in TS:
        if Sval(sts[t]) - S0 != t * e1:
            Sbad.append(u)
            break

say('knobs probed: %d   sample points t = %s' % (len(KNOB), TS))
say('S is exactly affine in every knob: %s   (exceptions: %s)' % (not Sbad, Sbad))
say('\nknobs with a non-affine check (these equations were DROPPED from the model):')
tot = set()
for u, s in sorted(bad_nonaffine.items()):
    say('   x_%-6d non-affine on %d checks %s' % (u, len(s), sorted(s)[:8]))
    tot |= s
say('union of non-affine checks: %d  %s' % (len(tot), sorted(tot)))

# what the earlier run had detected with only the 2-point probe
prev = set()
for u in KNOB:
    o = fv0.get(u, 0)
    s1 = st0.clone().set_free({u: o + 1})
    s2 = st0.clone().set_free({u: o + 2})
    for a in fr.checks:
        d = s1.av[a] - base[a]
        if s2.av[a] - base[a] != 2 * d:
            prev.add(a)
say('\ndetected by the original 2-point probe: %d  %s' % (len(prev), sorted(prev)))
missed = tot - prev
say('MISSED by the 2-point probe (would invalidate the model): %s' % sorted(missed))
if not missed:
    say('=> the two-point probe was sufficient; every kept row is exactly affine.')
    say('   The negative result stands as an exact statement over the knob class.')
else:
    say('=> the model was NOT exact; the negative result must be re-run excluding these.')
json.dump({'nonaffine_checks': sorted(tot), 'two_point_detected': sorted(prev),
           'missed': sorted(missed), 'S_affine': not Sbad},
          open(OD + '/fb_lincheck.json', 'w'))
say('DONE')
