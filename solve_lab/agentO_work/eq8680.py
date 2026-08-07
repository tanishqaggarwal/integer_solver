"""Is S = 0 FORCED, or merely satisfied at the witness?

eq8680 is the only equation containing a37887 = S^2.  If eq8680 has no other atom that anything
can move, then "eq8680 holds" implies S^2 = 0 implies S = 0 for EVERY satisfying assignment --
an unconditional lemma, not a knob-set-scoped one.  If it has other movable atoms, S != 0 is
compensable in principle and eq8680 stops being a 1-for-1 tax.

Dump the equation exactly, then test compensability.
"""
import sys, os, json, re, itertools
HERE = '/home/user/integer_solver/solve_lab/agentH_work'
sys.path.insert(0, HERE)
os.chdir(HERE)
sys.set_int_max_str_digits(20_000_000)
OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/eq8680.log', 'w', buffering=1)


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
assert st0.score() == 39026

E = 8680
m, sq, tl = frameB.eq_terms[E]
say('eq%d:  multiplier m = %s   squared-outer flag sq = %s   number of terms = %d'
    % (E, m, sq, len(tl)))
say('terms (coefficient, atom):')
for c, a in tl:
    src = frameB.atom_src[a] if a >= 0 else '<literal>'
    isck = a in set(fr.checks)
    val = st0.av.get(a, 0) if isck else None
    say('   c=%-6s atom %-6s  check=%-5s value=%s' % (c, a, isck, ('%d bits' % abs(val).bit_length()) if val else val))
    say('        src: %s' % src[:150])

A = 37887
say('\na%d is in equations: %s' % (A, fr.eq_of[A]))
say('current eq%d value: %s' % (E, st0.eq[E]))

# which atoms of eq8680 are checks (movable at all) vs definitions (identically satisfied)
ck = set(fr.checks)
inchk = [(c, a) for c, a in tl if a in ck]
notchk = [(c, a) for c, a in tl if a not in ck]
say('\nterms whose atom is a CHECK (can be nonzero): %s' % [(c, a) for c, a in inchk])
say('terms whose atom is a DEFINITION (identically zero in this frame): %d' % len(notchk))

say('\n--- VERDICT on whether S = 0 is forced')
if len(inchk) == 1 and inchk[0][1] == A:
    c = inchk[0][0]
    say('  eq%d reduces to  %s * %s * a%d = 0  with a%d = S^2.' % (E, m, c, A, A))
    say('  Every other atom in the equation is a DEFINITION and is identically zero in any')
    say('  assignment produced by this frame, so eq%d holds  <=>  S^2 = 0  <=>  S = 0.' % E)
    say('  ==> S = 0 is FORCED for every assignment that satisfies eq%d.' % E)
    say('      This is unconditional -- it does not depend on any knob set.')
else:
    say('  eq%d has %d check atoms; S != 0 is compensable in principle.' % (E, len(inchk)))
    say('  check atoms other than a%d: %s' % (A, [a for c, a in inchk if a != A]))

# how many free inputs can move each of the other check atoms of eq8680
say('\n--- can anything move the OTHER check atoms of eq%d?' % E)
for c, a in inchk:
    if a == A:
        continue
    sup = fr.SUPV.get(a, [])
    say('   a%-6d supported by %d free inputs; current value %s'
        % (a, len(sup), st0.av.get(a)))

json.dump({'m': str(m), 'sq': bool(sq), 'nterms': len(tl),
           'check_terms': [[str(c), a] for c, a in inchk]},
          open(OD + '/eq8680.json', 'w'))
say('DONE')
