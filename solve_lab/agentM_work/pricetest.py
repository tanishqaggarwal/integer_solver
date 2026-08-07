"""Validate the candidate-agnostic pricer on the known answer, and measure throughput."""
import sys, os, time, json, random
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import price as PR
import engine2 as E2

vd = PR.load_deliverable()
P = PR.Pricer(vd)
H4 = [642, 28730, 29854, 31864]          # P's four corrupted handles -- handles ONLY

print('=== VALIDATION: reproduce the deliverable from its four handles ===')
freed, demote = PR.closure(H4)
print(f'input handles      : {H4}')
print(f'closure -> freed   : {freed}')
print(f'closure -> demote  : {demote}')
print(f'engine2 PIN was    : {sorted(E2.PIN)}')
print(f'engine2 DEMOTE was : {sorted(E2.DEMOTE_ATOMS)}')
print(f'MATCH freed==PIN   : {sorted(freed) == sorted(E2.PIN)}')
print(f'MATCH demote       : {sorted(demote) == sorted(E2.DEMOTE_ATOMS)}')
print('  (the 5th demotion is DERIVED from the 4 handles, not supplied)')

vals = {u: vd[u] for u in freed}
t0 = time.time()
r = P.price_given(H4, vals)
t_given = time.time() - t0
print(f'\nprice_given -> score {r["score"]}, nbad {r["nbad"]}, {t_given:.2f}s')
print(f'  fails: {r["fails"]}')
print(f'  bad  : {r["bad"]}')
exp_fails = [12231, 12270, 12350, 14584, 18673, 22044, 29125]
exp_bad = [23616, 23617, 36659, 36660, 36661, 36662, 36663, 36664]
ok = (r['score'] == 39026 and r['fails'] == exp_fails and r['bad'] == exp_bad)
nd = sum(1 for u in range(PR.NV) if r['v'][u] != vd[u])
print(f'  vars differing from deliverable: {nd}')
print(f'\nVALIDATION {"PASSED" if ok and nd == 0 else "FAILED"}')

print('\n=== THROUGHPUT ===')
# path 1: values supplied
t0 = time.time(); n = 0
for _ in range(3):
    P.price_given(H4, vals); n += 1
t_given = (time.time() - t0) / n
print(f'  price_given  : {t_given:.2f}s per candidate  -> {3600/t_given:,.0f}/hour')

# path 2: values searched (the realistic path for a supplied handle set)
t0 = time.time()
rs = P.price_search(H4, maxk=2, budget=120)
t_search = time.time() - t0
print(f'  price_search : {t_search:.1f}s for the deliverable\'s own four')
print(f'     -> base(uncorrupted) {rs.get("base_score")}, best {rs.get("score")} via {rs.get("via")}, '
      f'{rs.get("nknobs")} affine knobs')

# random handle sets of the same size, to time the general case
prod = []
import re
for u in H.SEQ:
    a = H.definer[u][0]
    t = H.atoms[a]
    if re.fullmatch(r'x_%d - x_\d+ \* x_\d+' % u, t):
        prod.append(u)
print(f'\n  product-defined variables found in my frame: {len(prod)}')
rnd = random.Random(5)
times = []
for i in range(4):
    hs = rnd.sample(prod, 4)
    t0 = time.time()
    rr = P.price_search(hs, maxk=2, budget=60)
    dt = time.time() - t0
    times.append(dt)
    print(f'  random set {hs}: freed {len(rr.get("freed",[]) or [])}, '
          f'base {rr.get("base_score")}, best {rr.get("score")}, {dt:.1f}s')
if times:
    avg = sum(times) / len(times)
    print(f'\n  price_search average {avg:.1f}s per candidate -> ~{3600/avg:,.0f}/hour single-core')
