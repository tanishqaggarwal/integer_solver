"""Price O's delta0 in my frame, direct and with MINIMAL REPRESENTATIVES.

O's reductions, re-derived here so the compensations are explicit:

  atom 23616 = (x_7068 - x_2099) - 7376877*x_642
      shifting the external part by d and x_642 by (d_red - d)/7376877 leaves the atom
      unchanged, so d only matters mod 7376877 -> a 23-bit condition, not 2440-bit.
  atom 23618 = (x_4432 - x_19964) - x_28730
      external shift d, compensate x_28730 by -(d - d_red); x_28730 also enters
      x_28730 - p*x_9413 with +1, compensated by x_9413 += (d - d_red)/p.  So mod p.
  atom 36660 = 5113045*(x_7075*x_9118) - x_29854 : external moves by multiples of p,
      compensated in x_29854.  atom 36662 = x_7075*x_8731 likewise.

Why the reduction should matter: x_7068 touches only atoms 23616 and 34120, and atom 34120
drives 12 of the 25 baseline equations (check-in 42).  A 2440-bit move of x_7068 wrecks
34120; a 23-bit move is a far smaller perturbation of it.  That is the whole bet.
"""
import sys, os, json, time, collections
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine as EB, engine3 as E3
import price as PR, fscore

D = json.load(open('DELTA0_FOR_M.json'))
p = int(D['p'])
REG = D['region_equations']
BS = D['boundary_shifts']
Z = {int(k): int(v) for k, v in D['private_solution_z'].items()}
D4 = [642, 28730, 29854, 31864]

vd = PR.load_deliverable()
freed, demote = PR.closure(D4)
eng = E3.Eng(demote)
seed0 = {f: vd[f] for f in eng.FREE if vd[f] != 0}
v_base = eng.forward(seed0)
bad_base = eng.badatoms(v_base)
print(f'BASE = deliverable: score {fscore.score(bad_base)}  fails {sorted(fscore.fails(bad_base))}',
      flush=True)

FS = set(eng.FREE)
print(f'carriers free in my frame: '
      f'{ {u: (u in FS) for u in (7068, 4432, 9118, 8731, 642, 28730, 29854, 31864, 9413, 17325)} }',
      flush=True)

d1 = int(BS['23616']['carrier_increment'])      # via x_7068
d2 = int(BS['23618']['carrier_increment'])      # via x_4432
d3 = int(BS['36660']['carrier_increment'])      # via x_9118 (already /5113045)
d4 = int(BS['36662']['carrier_increment'])      # via x_8731
print(f'\ncarrier increments (bits): d1 {d1.bit_length()} d2 {d2.bit_length()} '
      f'd3 {d3.bit_length()} d4 {d4.bit_length()}', flush=True)


def bal(x, m):
    """balanced representative of x mod m, smallest magnitude"""
    r = x % m
    return r - m if r > m // 2 else r


def score_of(changes, label, dump=False):
    s = dict(seed0)
    for k, v in changes.items():
        if v:
            s[k] = v
        else:
            s.pop(k, None)
    v = eng.forward(s)
    av = eng.badatoms(v)
    fl = sorted(fscore.fails(av))
    sc = PR.NEQ - len(fl)
    reg_bad = [e for e in REG if e in fl]
    print(f'  [{label}] score {sc}  nbad {len(av)}  nfail {len(fl)}  '
          f'region still failing {len(reg_bad)}/13 {reg_bad}', flush=True)
    if dump or sc > 39026:
        fn = f'M_delta0_{label}_{sc}.json'
        json.dump({f"x_{i}": int(v[i]) for i in range(PR.NV) if v[i] != 0}, open(fn, 'w'))
        print(f'      wrote {fn}', flush=True)
    return sc, fl


print('\n=== A: delta0 applied directly (carriers shifted, z as INCREMENTS) ===', flush=True)
chA = {7068: vd[7068] + d1, 4432: vd[4432] + d2,
       9118: vd[9118] + d3, 8731: vd[8731] + d4}
for u, z in Z.items():
    if z and u in FS:
        chA[u] = vd[u] + z
score_of(chA, 'A_incr')

print('\n=== B: same carriers, z as ABSOLUTE values ===', flush=True)
chB = {7068: vd[7068] + d1, 4432: vd[4432] + d2,
       9118: vd[9118] + d3, 8731: vd[8731] + d4}
for u, z in Z.items():
    if u in FS:
        chB[u] = z
score_of(chB, 'B_abs')

print('\n=== C: MINIMAL REPRESENTATIVES ===', flush=True)
M1 = 7376877
d1r = bal(d1, M1)
k1 = (d1r - d1) // M1                      # x_642 compensation
d2r = bal(d2, p); k2 = d2 - d2r            # x_28730 compensation, multiple of p
d3r = bal(d3, p)
d4r = bal(d4, p)
print(f'  reduced (bits): d1 {d1.bit_length()}->{d1r.bit_length()} (mod {M1}), '
      f'd2 {d2.bit_length()}->{d2r.bit_length()}, d3 {d3.bit_length()}->{d3r.bit_length()}, '
      f'd4 {d4.bit_length()}->{d4r.bit_length()}', flush=True)
print(f'  d1 reduced value = {d1r}', flush=True)

for tag, use_z in (('C_incr', 'incr'), ('C_abs', 'abs')):
    ch = {7068: vd[7068] + d1r, 4432: vd[4432] + d2r,
          9118: vd[9118] + d3r, 8731: vd[8731] + d4r}
    for u, z in Z.items():
        if u in FS:
            ch[u] = (vd[u] + z) if use_z == 'incr' else z
    # compensations that keep the region atoms fixed under the reduction
    ch[642] = ch.get(642, vd[642]) + k1
    ch[28730] = ch.get(28730, vd[28730]) - k2
    if 9413 in FS and k2 % p == 0:
        ch[9413] = vd[9413] + (k2 // p)
    score_of(ch, tag)

print('\n=== D: reduce ONLY the x_7068 carrier (the 34120 lever), keep the rest as O gave ===',
      flush=True)
chD = {7068: vd[7068] + d1r, 4432: vd[4432] + d2,
       9118: vd[9118] + d3, 8731: vd[8731] + d4}
for u, z in Z.items():
    if z and u in FS:
        chD[u] = vd[u] + z
chD[642] = chD.get(642, vd[642]) + k1
score_of(chD, 'D_only7068')

print('\n=== E: carriers only, no private z (isolates the carrier cost) ===', flush=True)
score_of({7068: vd[7068] + d1r, 4432: vd[4432] + d2r,
          9118: vd[9118] + d3r, 8731: vd[8731] + d4r}, 'E_carriers_reduced')
score_of({9118: vd[9118] + d3r, 8731: vd[8731] + d4r}, 'E_freecarriers_only')
print('\nbaseline 39026 ; perfect 39033', flush=True)
