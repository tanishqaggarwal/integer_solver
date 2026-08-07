"""Fix coordinates and dimension before re-solving the lattice.

(1) Verify T's third calibration point in my scorer: deliverable with the 12 cofactors
    zeroed -> 39,021 / 12 failing, with T's exact list.  Price against 7 -> 12, not 7 -> 13.
(2) Establish the TRUE dimension of the space I am solving over.  Two classes:
      cofactors u        -- free in the base harness
      broken-atom wires  -- x642/x28730/x29854/x31864, NOT free in the base harness but
                            assignable here because the deliverable breaks their defining
                            atoms (which is exactly what engine3's demotion encodes)
    A coordinate is LIVE only if moving it actually changes the score/atoms.  Measured, not
    assumed, because eight of the twelve cofactors are already zero and do nothing.
"""
import sys, os, json, collections
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine as EB, engine3 as E3
import price as PR, fscore

COF12 = [105, 1329, 3387, 5081, 5676, 9413, 10903, 11436, 14393, 14768, 17325, 22820]
WIRES4 = [642, 28730, 29854, 31864]
CARRIERS = [7068, 4432, 9118, 8731]
T_EXPECT = [2554, 6816, 8124, 9123, 9421, 12231, 12270, 12350, 14584, 18673, 22044, 29125]

vd = PR.load_deliverable()
freed, demote = PR.closure(WIRES4)
eng = E3.Eng(demote)
seed = {f: vd[f] for f in eng.FREE if vd[f] != 0}
v0 = eng.forward(seed)
bad0 = eng.badatoms(v0)
F0 = sorted(fscore.fails(bad0))
print(f'deliverable in my frame: score {fscore.score(bad0)}, failing {F0}', flush=True)


def sc(changes):
    s = dict(seed)
    for k, v in changes.items():
        if v:
            s[k] = v
        else:
            s.pop(k, None)
    v = eng.forward(s)
    av = eng.badatoms(v)
    fl = sorted(fscore.fails(av))
    return PR.NEQ - len(fl), fl, av


print('\n=== (1) T calibration: deliverable with the 12 cofactors zeroed ===', flush=True)
s12, f12, a12 = sc({u: 0 for u in COF12})
print(f'  score {s12}   failing {len(f12)}: {f12}', flush=True)
print(f'  T expects 39021 / 12 / {T_EXPECT}', flush=True)
print(f'  MATCH score {s12 == 39021}   MATCH list {f12 == T_EXPECT}', flush=True)
print(f'  support (bad atoms) {sorted(a12)}', flush=True)
print(f'  deliverable support {sorted(bad0)}   support identical: '
      f'{sorted(a12) == sorted(bad0)}', flush=True)
print(f'\n  ==> PRICE AGAINST 7 -> {len(f12)}  (gap {len(f12) - len(F0)})', flush=True)

print('\n=== (2a) which of the 12 cofactors are LIVE (zero each alone) ===', flush=True)
live_cof = []
for u in COF12:
    s1, f1, _ = sc({u: 0})
    d = len(f1) - len(F0)
    tag = 'LIVE' if d else 'dead (already 0)' if vd[u] == 0 else 'dead'
    if d:
        live_cof.append(u)
    print(f'  x_{u:<6d} value {"0" if vd[u]==0 else str(vd[u].bit_length())+" bits":>10s}  '
          f'zeroing -> failures {len(f1):3d}  delta +{d}   {tag}', flush=True)
print(f'  LIVE cofactors: {live_cof}  -> dimension {len(live_cof)}', flush=True)

print('\n=== (2b) the four broken-atom wires: free here? and live? ===', flush=True)
FSbase = set(EB.FREE)
FShere = set(eng.FREE)
live_wire = []
for u in WIRES4:
    s1, f1, _ = sc({u: 0})
    d = len(f1) - len(F0)
    if d:
        live_wire.append(u)
    print(f'  x_{u:<6d} free in BASE harness: {u in FSbase:5}  free in MY frame: '
          f'{u in FShere:5}  value {vd[u].bit_length():4d} bits  zeroing -> +{d}', flush=True)
print(f'  LIVE wires: {live_wire}  -> dimension {len(live_wire)}', flush=True)

print('\n=== (2c) carriers ===', flush=True)
live_car = []
for u in CARRIERS:
    s1, f1, _ = sc({u: vd[u] + 1})
    d = len(f1) - len(F0)
    if d:
        live_car.append(u)
    print(f'  x_{u:<6d} free in BASE: {u in FSbase:5}  free here: {u in FShere:5}  '
          f'+1 -> failures {len(f1):3d} (delta {d:+d})', flush=True)

print('\n=== TRUE DIMENSION ===', flush=True)
print(f'  live cofactors {len(live_cof)}: {live_cof}')
print(f'  live wires     {len(live_wire)}: {live_wire}')
print(f'  carriers        {len(CARRIERS)}: {CARRIERS}')
tot = sorted(set(live_cof) | set(live_wire) | set(CARRIERS))
print(f'  UNION: {len(tot)} coordinates -> {tot}')
json.dump({'live_cofactors': live_cof, 'live_wires': live_wire,
           'carriers': CARRIERS, 'coords': tot,
           'far_side_failures': f12, 'far_side_score': s12},
          open('dimcheck.json', 'w'), indent=1)
print('\nwrote dimcheck.json', flush=True)
