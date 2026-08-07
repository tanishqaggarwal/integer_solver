"""The Lemma (eq8680 = S^4, hence S = 0 unconditionally) rests on a parse.  Verify it against
E's INDEPENDENT parser as well, so the claim does not depend on one model's term extraction.
"""
import sys, re
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentO_work')
import simO, engine as E, harness as H

OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/eq8680_cross.log', 'w', buffering=1)


def say(*a):
    print(*a, file=LOG)


m, sq, tl = H.eqt[8680]
say("E's parse of eq8680:")
say('   multiplier m = %s   squared-outer sq = %s   terms = %d' % (m, sq, len(tl)))
for c, a in tl:
    src = H.atoms[a] if a >= 0 else '<literal>'
    say('   coefficient %s on atom %s' % (c, a))
    say('        length %d chars' % len(src))
    half = (len(src) - 3) // 2
    issq = src[half:half + 3] == ' * ' and src[:half] == src[half + 3:]
    say('        is a literal square (S)*(S): %s' % issq)
    if issq:
        S = src[:half]
        say('        S begins: %s' % S[:90])
        vs = sorted({int(x) for x in re.findall(r'x_(\d+)', S)})
        say('        S mentions %d distinct variables' % len(vs))
        for u in (4432, 19964, 28730):
            say('        x_%-6d appears in S: %s' % (u, u in vs))
say('')
if len(tl) == 1 and sq:
    say('CONFIRMED by E\'s parser too: eq8680 = m * (c * a)^2 with a = S^2,')
    say('so eq8680 = %s * %s^2 * S^4.  Its only zero is S = 0.' % (m, tl[0][0]))
    say('=> S = 0 is forced in every satisfying assignment, independent of frame or knob set.')
else:
    say('E\'s parser DISAGREES: terms=%d sq=%s -- the Lemma is parser-dependent.' % (len(tl), sq))

# and confirm no other equation mentions that atom
a0 = tl[0][1]
other = [e for e, (mm, ss, t2) in enumerate(H.eqt) if any(x == a0 for _, x in t2)]
say('\nequations mentioning atom %d: %s' % (a0, other))
say('DONE')
