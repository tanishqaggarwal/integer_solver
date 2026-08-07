"""Can re-orientation retire any of the 7 failing equations outright?

A frame's DEFINITION atoms are identically zero for every choice of free inputs (that is what
`x_t := rest` means).  Only CHECK atoms carry value into the equations.  So an equation all of
whose atoms are definitions is satisfied unconditionally, and re-orienting a check atom into a
definition retires every equation that atom alone was breaking -- at the price of turning whatever
used to define that variable into a new check.

This asks the sharp question: for each of the 7 equations failing at the 39,026 witness, which of
its atoms are checks, do those checks admit a legal unit target, and what does taking it cost?
"""
import ast, json, pickle, os
from collections import defaultdict
import model, ev
from orient import unit_targets

HERE = os.path.dirname(os.path.abspath(__file__))
d = model.get()
atom_src = d['atom_src']
atom_vars = d['atom_vars']
eq_terms = d['eq_terms']
NA = len(atom_src)
NV = 38748
F = pickle.load(open(os.path.join(HERE, 'fwd2.pkl'), 'rb'))
definer = F['definer']
checks = set(F['checks'])
free0 = set(F['free0'])
tgt = F['tgt']
defatom_of = {}                      # atom -> variable it defines
for v in range(NV):
    if definer[v] >= 0:
        defatom_of[definer[v]] = v

FAIL = [12231, 12270, 12350, 14584, 18673, 22044, 29125]

# witness atom values
import optN
from optN import make, POOL
stw = make(POOL)
print('witness score %d, failing %s' % (stw.score(), sorted(stw.fails)))

alts = {}
print('\n=== the 7 failing equations, atom by atom ===')
freeable = defaultdict(list)
for e in FAIL:
    m, sq, tl = eq_terms[e]
    print('\neq %-6d  m=%s sq=%s  %d atom terms' % (e, m, sq, len(tl)))
    for c, a in tl:
        isck = a in checks
        val = stw.av.get(a, 0)
        u = alts.setdefault(a, unit_targets(atom_src[a]))
        state = 'CHECK' if isck else 'def(x_%d)' % defatom_of.get(a, -1)
        cand = []
        for v, s in u.items():
            if v in free0:
                where = 'free input'
            elif definer[v] >= 0:
                where = 'defined by atom %d' % definer[v]
            else:
                where = 'unresolved'
            cand.append('x_%d(sign%+d, %s)' % (v, s, where))
        print('   coef %-4s atom %-6d %-14s witness value %s   unit targets: %s'
              % (c, a, state, ('0' if val == 0 else 'NONZERO(%d digits)' % len(str(abs(val))))),
              end='')
        print(' %s' % ('; '.join(cand) if cand else 'NONE'))
        if isck and val != 0 and u:
            freeable[e].append((a, u))

print('\n=== summary ===')
tot = 0
for e in FAIL:
    m, sq, tl = eq_terms[e]
    cks = [a for c, a in tl if a in checks]
    nz = [a for c, a in tl if a in checks and stw.av.get(a, 0) != 0]
    orientable = [a for a in nz if alts.get(a)]
    tot += len(orientable)
    print('eq %-6d: %d atoms, %d checks, %d nonzero at witness, %d of those re-orientable %s'
          % (e, len(tl), len(cks), len(nz), len(orientable), orientable))
print('\ntotal re-orientable nonzero check atoms across the 7 failing equations: %d' % tot)
