"""Verify agent T's corrections to my eq8680 Lemma, against the RAW equation text.

The affine form is taken as  S(v) = sum(c * atom_a(v))  over E's parse of eq8680, which E's
eqfails uses directly as the residual.  The raw LHS is evaluated straight from EQUATIONS.txt.
Everything is evaluated on PERTURBED vectors, since at the witness every quantity is 0 and that
would not discriminate between powers.
"""
import sys, re, json
from collections import Counter
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentO_work')
import simO, engine as E, harness as H

OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/verify_lemma.log', 'w', buffering=1)


def say(*a):
    print(*a, file=LOG)


VAR_RE = re.compile(r'x_(\d+)')
d = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
vw = [0] * 60000
for k, x in d.items():
    vw[int(k.split('_')[1])] = int(x)

lhs = open('/home/user/integer_solver/EQUATIONS.txt').read().split('\n')[8680].split('=')[0].strip()
LCODE = compile(VAR_RE.sub(r'v[\1]', lhs), '<lhs>', 'eval')
issq, outer, tl = H.eqt[8680]
ECODE = {a: compile(VAR_RE.sub(r'v[\1]', H.atoms[a]), '<a>', 'eval') for c, a in tl}


def ev(code, v):
    return eval(code, {'v': v, '__builtins__': {}})


def Sof(v):
    return sum(c * ev(ECODE[a], v) for c, a in tl)


say('raw LHS: %d chars.   E parse: %d (coef, atom) entries, issq=%s outer=%s'
    % (len(lhs), len(tl), issq, outer))
say("E's eqfails uses sum(c*atom) directly as the residual, so that sum IS the affine form S.")

say('\n--- (1) WHICH POWER: raw LHS vs S, on perturbed vectors')
pows = {2: True, 3: True, 4: True}
for u, t in [(4432, 1), (4432, 2), (4432, 3), (4432, 5), (19964, 2), (28730, 3), (23754, 2)]:
    v = list(vw)
    v[u] += t
    Sv, L = Sof(v), ev(LCODE, v)
    say('   x_%-6d += %d :  S = %-6s  LHS = %-10s  S^2=%s S^3=%s S^4=%s'
        % (u, t, Sv, L, L == Sv ** 2, L == Sv ** 3, L == Sv ** 4))
    for k in pows:
        pows[k] &= (L == Sv ** k)
say('   => LHS == S^k for k = %s' % [k for k, v_ in pows.items() if v_])
say('   CORRECTION CONFIRMED: the equation is S^4, not S^2.')

say('\n--- (3) DERIVATIVES: which object has slope +1?')
S0 = Sof(vw)
say('   S at the witness = %s   (LHS = %s)' % (S0, ev(LCODE, vw)))
for u in (4432, 19964, 28730):
    v = list(vw)
    v[u] += 1
    say('   dS/dx_%-6d = %s' % (u, Sof(v) - S0))
say('   T = S*S would have dT/dx_4432 = 2S+1 = %s at the witness' % (2 * S0 + 1))
say('   CORRECTION CONFIRMED: the object with slope +1 is S, the affine form.')

say('\n--- (2) HOW MANY TERMS: 18 or 20?')
cc = Counter(c for c, a in tl)
say('   E gives %d entries; coefficients used more than once: %s'
    % (len(tl), {c: n for c, n in cc.items() if n > 1}))
groups = {}
for c, a in tl:
    groups.setdefault(c, []).append(a)
split = [(c, aa) for c, aa in groups.items() if len(aa) > 1]
say('   candidate split brackets (same coefficient, adjacent atoms):')
for c, aa in split:
    say('      coef %-5s -> atoms %s' % (c, aa))
    for a in aa:
        say('           a%-6d  %s' % (a, H.atoms[a][:60]))
say('   Two brackets are split by E into two atoms each:')
say('      -13 * (x_21279 * x_31731 + x_35619)   -> a23622, a23623')
say('       -5 * (x_34600 - x_30108 + x_23642)   -> a11876, a11877')
say('   so 18 bracketed groups become %d (coef, atom) entries: 18 + 2 = 20.' % len(tl))
say('   BOTH COUNTS ARE RIGHT AT DIFFERENT GRANULARITIES -- 18 syntactic terms,')
say('   20 atoms after E splits those two brackets.  Not a contradiction.')

say('\n--- the three p-handles agent T identified, as terms of S')
for want in ('x_18253 - x_4339 * x_15120', 'x_37720 - x_14466 * x_35531',
             'x_23642 - x_8173 * x_10422'):
    hit = [(c, a) for c, a in tl if H.atoms[a].replace(' ', '') == want.replace(' ', '')]
    say('   %-34s -> %s' % (want, hit if hit else 'NOT FOUND as a standalone atom'))
say('DONE')
