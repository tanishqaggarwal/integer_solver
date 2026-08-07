"""Re-run the incidence filter over the FULL p-handle family, not the 3,681 census.

T's correction: a p-handle is defined by `h = p*u` with `u` free.  L's census was
delimited by GUARD SHAPE (slot links), so it omitted p-handles whose guards are stage
checks or leaf pins -- 33 of them, three of which are incident.  The defining property is
the algebraic form, not the guard, so I enumerate on the form directly:

    definer atom of the shape  x_h - x_i * x_j   with  value(x_i) == p  and  x_j free
                                                  (or the mirrored operand order)

and then test incidence EXACTLY -- atom a is in equation e iff a appears in eqt[e]'s terms
-- rather than via the cofactor-marker shortcut, which is what let the earlier pass miss
bare- and linearly-defined handles.
"""
import sys, os, re, json, collections
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine as EB
import price as PR, fscore

P = 115792089237316195423570985008687907853269984665640564039457584007908834671663
VD = PR.load_deliverable()
FREE = set(EB.FREE)

D7 = [12231, 12270, 12350, 14584, 18673, 22044, 29125]
T12 = [2554, 6816, 8124, 9123, 9421, 12231, 12270, 12350, 14584, 18673, 22044, 29125]
M25 = json.load(open('baseline_sets.json'))['A']

# ---------- the full p-handle family, by algebraic form ----------
PROD = re.compile(r'x_(\d+) - x_(\d+) \* x_(\d+)$')
family = {}          # atom -> (h, u, other)
both = 0
# scan ALL atoms, not just definer atoms: a p-handle atom need not be the definer of its
# own h (x_23642's definer is the BARE atom `x_23642`, while its p-handle atom is a
# separate check atom).  Restricting to definers is what undercounted this to 1,256.
for a in range(len(H.atoms)):
    m = PROD.fullmatch(H.atoms[a])
    if not m:
        continue
    h, i, j = int(m.group(1)), int(m.group(2)), int(m.group(3))
    opts = []
    if VD[i] == P and j in FREE:
        opts.append((j, i))
    if VD[j] == P and i in FREE:
        opts.append((i, j))
    if len(opts) == 2:
        both += 1
    if opts:
        family[a] = (h, opts[0][0], opts[0][1])
print(f'p-handle family by algebraic form: {len(family)} atoms '
      f'({len(family) + both} counting both operand orders)', flush=True)
print(f'  T reports 3,707 (3,714 both orders)', flush=True)

# ---------- exact incidence ----------
def incident(targets, label):
    hit = {}
    for a, (h, u, pv) in family.items():
        es = [e for e in targets
              if any(aa == a for c, aa in H.eqt[e][2] if aa >= 0)]
        if es:
            hit[a] = (h, u, es)
    print(f'\n=== incident against {label} ({len(targets)} equations): '
          f'{len(hit)} atoms ===', flush=True)
    for a, (h, u, es) in sorted(hit.items(), key=lambda kv: -len(kv[1][2])):
        print(f'  a{a:<6d} h=x{h:<6d} u=x{u:<6d}  rt {len(es):2d}  {es}', flush=True)
    return hit


h7 = incident(D7, "the deliverable's 7 failing")
h12 = incident(T12, "T's 12-equation far side")
h25 = incident(M25, 'my 25-equation uncorrupted baseline')

# ---------- cross-check T's three ----------
TSU = [10422, 15120, 35531]
print('\n=== cross-check: T\'s three newly-found incident handles ===', flush=True)
byu = {u: (a, h) for a, (h, u, _) in family.items()}
for u in TSU:
    infam = u in byu
    a = byu.get(u, (None, None))[0]
    in7 = a in h7 if a else False
    in12 = a in h12 if a else False
    print(f'  u=x{u:<6d} in family: {infam}   atom a{a}   incident to the 7: {in7}   '
          f'to the 12: {in12}', flush=True)

# ---------- verify the criterion T checked: eqs(u) == eqs(atom_u) ----------
print('\n=== verify T\'s criterion eqs(u) == eqs(atom_u) over the whole family ===',
      flush=True)
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
import checker
codes, varsets = checker.load_equations()
atom_eqs = collections.defaultdict(set)
for e, (issq, outer, terms) in enumerate(H.eqt):
    for c, a in terms:
        if a >= 0:
            atom_eqs[a].add(e)
var_eqs = collections.defaultdict(set)
for e, vs in enumerate(varsets):
    for v in vs:
        var_eqs[v].add(e)
viol = 0; checked = 0
for a, (h, u, pv) in family.items():
    checked += 1
    if var_eqs.get(u, set()) != atom_eqs.get(a, set()):
        viol += 1
print(f'  checked {checked}, violations {viol}  '
      f'({"criterion holds" if viol == 0 else "CRITERION FAILS"})', flush=True)

json.dump({'family_size': len(family), 'both_orders': len(family) + both,
           'incident_7': {str(a): {'h': v[0], 'u': v[1], 'eqs': v[2]} for a, v in h7.items()},
           'incident_12': {str(a): {'h': v[0], 'u': v[1], 'eqs': v[2]} for a, v in h12.items()},
           'incident_25': {str(a): {'h': v[0], 'u': v[1], 'eqs': v[2]} for a, v in h25.items()},
           'criterion_violations': viol},
          open('pfamily.json', 'w'), indent=1)
print(f'\nwrote pfamily.json', flush=True)
print(f'\nSPACE: 2^{len(h12)} = {2**len(h12):,} (against the 12-equation far side)', flush=True)
