"""Cross-check against agent O's Lemma: is T = 0 what my 924/924 obstruction measures?

O: eq8680 = T^2 with T a linear form, dT/dx_4432 = +1, dT/dx_28730 = -1, zero elsewhere; a square
has a single zero locus so every satisfying assignment has T = 0.

In my frame eq_terms[8680] = (m=1, sq=True, [(1, 37887)]), so the inner linear form is exactly the
check atom 37887 -- and `optN.inner` returns the INNER form, never its square, so my linear model
already carries T rather than T^2.  atom 37887 is one of pool.py's REGION_ATOMS.

Question: does T = 0 explain the 924/924 p-obstruction, or is it a separate constraint?
"""
import json, itertools
import ev, model, optN, zsolve
from optN import make, build, inner, WIT, POOL

d = model.get()
print('=== O\'s Lemma in this frame ===')
m, sq, tl = ev.eq_terms[8680]
print('eq_terms[8680] = (m=%s, sq=%s, terms=%s)' % (m, sq, tl))
src = d['atom_src'][37887]
print('atom 37887 source (T), %d chars:' % len(src))
print('   ', src[:400] + ('...' if len(src) > 400 else ''))
av = d['atom_vars'][37887]
print('T depends on %d variables; 4432 in T: %s ; 28730 in T: %s' % (len(av), 4432 in av, 28730 in av))

# numeric derivative of T w.r.t. each of its variables, from a base state
st0 = make([])
base = st0.av.get(37887, 0)
print('\nT at make([])       = %s' % base)
derivs = {}
for v in sorted(av):
    h = st0.clone()
    try:
        h.set_free({v: st0.v[v] + 1})
    except Exception:
        continue
    dd = h.av.get(37887, 0) - base
    if dd:
        derivs[v] = dd
print('numeric dT/dx over the variables that move it: %s' % derivs)

stw = make(POOL)
print('T at make(POOL) (=witness) = %s' % stw.av.get(37887, 0))
st28 = make([28730])
print('T at make([28730])         = %s' % st28.av.get(37887, 0))

print('\n=== which equation does detaching 28730 actually buy? ===')
f0 = set(st0.fails)
f1 = set(st28.fails)
print('make([]) failing %d, make([28730]) failing %d' % (len(f0), len(f1)))
print('fixed by detaching 28730: %s' % sorted(f0 - f1))
print('broken by detaching 28730: %s' % sorted(f1 - f0))

print('\n=== is 8680 in the regions? ===')
for D, tag in (([], 'make([])'), ([28730], 'make([28730])'), (WIT, 'witness')):
    b = build(make(list(D)))
    print('%-16s |R|=%d   8680 in R: %-5s   T=inner(8680)=%s'
          % (tag, len(b['R']), 8680 in b['R'], inner(make(list(D)), 8680)))

print('\n=== does T=0 explain the 924/924 obstruction? ===')
print('The witness region is the 12 rows EXCLUDING 8680, so T=0 already holds there and is not')
print('among the constraints the 924 six-row subsets are being asked to satisfy.  Test directly:')
st = make(WIT)
b = build(st)
R, M, bb, n = b['R'], b['M'], b['b'], b['n']
print('witness region R = %s' % (R,))
print('8680 present: %s' % (8680 in R))

# Now the 13-row region: is 8680 zeroable there, and is its obstruction also p?
st13 = make([642])
b13 = build(st13)
R13 = b13['R']
i8680 = R13.index(8680)
Z = zsolve.ZSolver(b13['M'], b13['b'], b13['n'])
print('\n13-row region (D=[642]): 8680 is row %d, current value T=%s' % (i8680, b13['b'][i8680]))
print('  row 8680 individually integrally zeroable? %s' % Z.solvable([i8680]))
opt_all, rows_all, _, _ = zsolve.max_zero_rows(b13['M'], b13['b'], b13['n'], len(R13))
print('  OPT=%d, an optimal row set = %s (equations %s)'
      % (opt_all, rows_all, [R13[i] for i in rows_all]))
print('  is 8680 in some optimal set? %s' % (i8680 in rows_all))
# best achievable when 8680 is forced zero
best = 0
for size in range(len(R13), 0, -1):
    found = False
    for A in itertools.combinations(range(len(R13)), size):
        if i8680 in A and Z.solvable(A):
            best = size
            found = True
            break
    if found:
        break
print('  max zeroable rows among the 13 SUBJECT TO 8680 being zeroed: %d' % best)
