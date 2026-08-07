"""Verify the two congruences have EXACTLY the claimed shape at all 383 stages, structurally,
from the raw atom definitions -- no trust in anyone's prose."""
import sys, os, json, re, ast
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
import model
d = model.get(); atom_src = d['atom_src']; atom_vars = d['atom_vars']
V = re.compile(r'x_(\d+)')
stages = json.load(open('/home/user/integer_solver/solve_lab/agentQ_work/qstages.json'))['stages']
print('stages:', len(stages))

# index: for each variable, the atoms of the form  x_t - RHS  that DEFINE it
defs = {}
for a, s in enumerate(atom_src):
    t = ast.parse(s, mode='eval').body
    if isinstance(t, ast.BinOp) and isinstance(t.op, ast.Sub) and isinstance(t.left, ast.Name):
        defs.setdefault(int(t.left.id[2:]), []).append((a, ast.unparse(t.right)))

def uniq_def(v):
    L = defs.get(v, [])
    return L[0][1] if len(L) == 1 else None

def N(v): return 'x_%d' % v

ok1 = ok2 = 0; K_wires = set(); bad = []
for i, g in enumerate(stages):
    ua, ub, ya, yb, u3, y3 = g['ua'], g['ub'], g['ya'], g['yb'], g['u3'], g['y3']
    dx, dy, S, R1 = g['dx'], g['dy'], g['S'], g['R1']
    # --- congruence 1 ---
    e = {}
    e['dx'] = uniq_def(dx); e['dy'] = uniq_def(dy); e['S'] = uniq_def(S); e['R1'] = uniq_def(R1)
    good = (e['dx'] == '%s - %s' % (N(ua), N(ub)) and e['dy'] == '%s - %s' % (N(ya), N(yb)))
    # R1 = (S*dx^2) - (dy^2) : resolve the two operands
    m = re.fullmatch(r'x_(\d+) - x_(\d+)', e['R1'] or '')
    if good and m:
        p1, p2 = int(m.group(1)), int(m.group(2))
        d1, d2 = uniq_def(p1), uniq_def(p2)
        m1 = re.fullmatch(r'x_(\d+) \* x_(\d+)', d1 or '')
        if m1 and d2 == '%s * %s' % (N(dy), N(dy)):
            f1, f2 = int(m1.group(1)), int(m1.group(2))
            sq = uniq_def(f2) if f1 == S else (uniq_def(f1) if f2 == S else None)
            if sq == '%s * %s' % (N(dx), N(dx)):
                # S = ua + ub + u3 + K  (chase the two-term adds)
                chain, cur, guard = [], e['S'], 0
                while cur and guard < 6:
                    mm = re.fullmatch(r'x_(\d+) \+ x_(\d+)', cur)
                    if not mm: break
                    l, r = int(mm.group(1)), int(mm.group(2))
                    chain.append(r); cur = uniq_def(l); guard += 1
                    if cur is None: chain.append(l); break
                terms = set(chain)
                if {ua, ub, u3} <= terms and len(terms) == 4:
                    K_wires |= (terms - {ua, ub, u3}); ok1 += 1
                else: bad.append((i, 'S terms', sorted(terms)))
            else: bad.append((i, 'dx^2', sq))
        else: bad.append((i, 'R1 operands', d1, d2))
    else: bad.append((i, 'R1/dx/dy', e['R1'], e['dx'], e['dy']))
    # --- congruence 2: find the wire  A*dx  and  B*dy , then their difference ---
    Adx = [v for v, L in defs.items() for _, r in L
           if r in ('x_%d * x_%d' % (0, 0),)]  # placeholder, real search below
    cand = []
    for v, L in defs.items():
        for _, r in L:
            mm = re.fullmatch(r'x_(\d+) \* x_(\d+)', r)
            if mm and {int(mm.group(1)), int(mm.group(2))} == {dx} | ({int(mm.group(1))} ^ {int(mm.group(2))}) - {dx}:
                pass
    # simpler: wires defined as (something * dx) and (dy * something)
    P = [v for v, L in defs.items() for _, r in L
         if re.fullmatch(r'x_(\d+) \* x_(\d+)', r) and dx in [int(z) for z in V.findall(r)]
         and len(set(V.findall(r))) == 2]
    Q = [v for v, L in defs.items() for _, r in L
         if re.fullmatch(r'x_(\d+) \* x_(\d+)', r) and dy in [int(z) for z in V.findall(r)]
         and len(set(V.findall(r))) == 2]
    hit = None
    for v, L in defs.items():
        for _, r in L:
            mm = re.fullmatch(r'x_(\d+) - x_(\d+)', r)
            if mm and int(mm.group(1)) in P and int(mm.group(2)) in Q:
                pa, pb = int(mm.group(1)), int(mm.group(2))
                ra, rb = uniq_def(pa), uniq_def(pb)
                Av = [int(z) for z in V.findall(ra) if int(z) != dx]
                Bv = [int(z) for z in V.findall(rb) if int(z) != dy]
                if Av and Bv:
                    da, db = uniq_def(Av[0]), uniq_def(Bv[0])
                    if da == '%s + %s' % (N(y3), N(yb)) and db == '%s - %s' % (N(ub), N(u3)):
                        hit = (v, pa, pb); break
        if hit: break
    if hit: ok2 += 1
    else: bad.append((i, 'R2 not matched'))

print('congruence 1  (S*dx^2 - dy^2, S = ua+ub+u3+K) matched at %d / %d stages' % (ok1, len(stages)))
print('congruence 2  ((y3+yb)*dx - (ub-u3)*dy)      matched at %d / %d stages' % (ok2, len(stages)))
print('distinct K wires across all stages: %d -> %s' % (len(K_wires), sorted(K_wires)[:6]))
for b in bad[:8]: print('  MISMATCH', b)
json.dump({'ok1': ok1, 'ok2': ok2, 'n': len(stages), 'K': sorted(K_wires)}, open('w_shape.json','w'))
