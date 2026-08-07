"""Verify agent T's three corrections to my eq8680 Lemma against the RAW equation text,
using no parser at all.  At the witness every quantity is 0, which would not discriminate,
so everything is evaluated on perturbed vectors.

Claims:
  (1) the equation is S^4, not S^2 (two nesting levels: LHS = T*T, T = S*S)
  (2) the affine form has 18 terms, not the 20 I reported
  (3) the object with derivative +1 is S; dT/dx_4432 = 2S+1, not 1
"""
import sys, re, json
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentO_work')
import simO, engine as E, harness as H

OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/verify_lemma.log', 'w', buffering=1)


def say(*a):
    print(*a, file=LOG)


VAR_RE = re.compile(r'x_(\d+)')
NV = 60000
d = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
vw = [0] * NV
for k, x in d.items():
    vw[int(k.split('_')[1])] = int(x)

line = open('/home/user/integer_solver/EQUATIONS.txt').read().split('\n')[8680]
lhs = line.split('=')[0].strip()
say('raw line 8680, LHS length %d chars' % len(lhs))


def peel(s):
    """If s is exactly (X) * (X) with identical factors, return X, else None."""
    s = s.strip()
    h = (len(s) - 3) // 2
    if s[h:h + 3] == ' * ' and s[:h] == s[h + 3:]:
        return s[:h]
    return None


lv = [lhs]
while True:
    p = peel(lv[-1])
    if p is None:
        break
    lv.append(p)
say('nesting depth: %d  (lengths %s)' % (len(lv) - 1, [len(x) for x in lv]))
S_txt = lv[-1]
say('innermost form S: %d chars' % len(S_txt))

CODE = {i: compile(VAR_RE.sub(r'v[\1]', t), '<x>', 'eval') for i, t in enumerate(lv)}


def ev(code, v):
    return eval(code, {'v': v, '__builtins__': {}})


say('\n--- (1) which power?  evaluate LHS and S on perturbed vectors')
ok4 = ok2 = True
for t in (1, 2, 3, 5):
    v = list(vw)
    v[4432] += t
    Sv = ev(CODE[len(lv) - 1], v)
    L = ev(CODE[0], v)
    say('   x_4432 += %d :  S = %s   LHS = S^4 ? %s   LHS = S^2 ? %s'
        % (t, Sv, L == Sv ** 4, L == Sv ** 2))
    ok4 &= (L == Sv ** 4)
    ok2 &= (L == Sv ** 2)
say('   => LHS == S^4 for all samples: %s      LHS == S^2: %s' % (ok4, ok2))

say('\n--- (3) derivatives: which object has slope +1?')
S0 = ev(CODE[len(lv) - 1], vw)
say('   S at the witness = %s' % S0)
for u in (4432, 19964, 28730):
    v = list(vw)
    v[u] += 1
    say('   dS/dx_%-6d = %s' % (u, ev(CODE[len(lv) - 1], v) - S0))
if len(lv) >= 3:
    T0 = ev(CODE[len(lv) - 2], vw)
    v = list(vw)
    v[4432] += 1
    say('   dT/dx_4432   = %s   (T = S*S, so this is 2S+1 = %s)'
        % (ev(CODE[len(lv) - 2], v) - T0, 2 * S0 + 1))

say('\n--- (2) how many terms?  18 or 20?')
issq, outer, tl = H.eqt[8680]
say("   E's parse gives %d (coef, atom) entries" % len(tl))
# top-level '+' split of S, respecting parentheses
parts, depth, cur = [], 0, ''
for ch in S_txt:
    if ch == '(':
        depth += 1
    elif ch == ')':
        depth -= 1
    if depth == 0 and ch == '+' :
        parts.append(cur)
        cur = ''
    else:
        cur += ch
parts.append(cur)
say('   raw text of S splits into %d top-level "+"-separated groups' % len(parts))
for p in parts[:22]:
    say('        %s' % p.strip()[:78])

# do E's 20 entries sum to S?
ECODE = {a: compile(VAR_RE.sub(r'v[\1]', H.atoms[a]), '<a>', 'eval') for c, a in tl}
say('\n   does sum(c * atom) over E\'s entries equal S ?')
allok = True
for t in (0, 1, 2, 5):
    v = list(vw)
    v[4432] += t
    tot = sum(c * ev(ECODE[a], v) for c, a in tl)
    Sv = ev(CODE[len(lv) - 1], v)
    say('      x_4432 += %d :  sum = %s   S = %s   equal: %s' % (t, tot, Sv, tot == Sv))
    allok &= (tot == Sv)
say('      => E\'s entries sum exactly to S: %s' % allok)

from collections import Counter
cc = Counter(c for c, a in tl)
dupes = {c: n for c, n in cc.items() if n > 1}
say('\n   coefficients appearing more than once in E\'s parse: %s' % dupes)
say('   (a bracket like  -13 * (A + B)  becomes TWO entries sharing coefficient -13,')
say('    which is how 18 bracketed groups become 20 (coef, atom) entries.)')
say('   %d groups + %d extra splits = %d entries'
    % (len(parts), len(tl) - len(parts), len(tl)))
say('DONE')
