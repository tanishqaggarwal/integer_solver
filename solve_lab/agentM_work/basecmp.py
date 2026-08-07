"""Is the 25-equation filter baseline-INDEPENDENT, or an artifact of E's orientation?

Baseline A (what I used): E's ORIGINAL orientation, full forward from the deliverable's
  free inputs. Every definer atom is forced to zero -- including the 5 the deliverable
  needs nonzero. This is the known-defective frame.

Baseline B (the deliverable's own): start from the deliverable's actual vector in the
  CORRECTED engine and UN-corrupt it in place -- set each freed variable to the value its
  own definer atom prescribes, evaluated at the deliverable's state -- then re-propagate.
  Cofactors and everything else keep the deliverable's values.

If A and B give the same failing-equation set, the filter is a property of the instance at
this free-input configuration. If they differ, it is a property of the orientation and L
must not discard candidates on it.
"""
import sys, os, collections, json, math
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine as EB, engine3 as E3
import price as PR, fscore

vd = PR.load_deliverable()
h0 = [642, 28730, 29854, 31864]
freed, demote = PR.closure(h0)
eng = E3.Eng(demote)
print(f'freed {freed}\ndemote {demote}\n')

# ---------- Baseline A : E's original orientation ----------
seed_unc = {f: vd[f] for f in EB.FREE if vd[f] != 0}
vA = EB.forward(seed_unc)
badA = EB.badatoms(vA)
FA = sorted(fscore.fails(badA))
print(f'BASELINE A (E orientation, full forward): score {fscore.score(badA)}, '
      f'{len(FA)} failures, {len(badA)} bad atoms')

# ---------- Baseline B : un-corrupt the deliverable in its own frame ----------
# For each freed var u with definer atom i, solve atom i == 0 for u AT THE DELIVERABLE'S
# STATE, i.e. put u back on the value its own definition prescribes. Iterate to a fixpoint
# because x_7068's definer references x_642.
v = list(vd)
ns = {'v': v, '__builtins__': {}}
order = sorted(freed, key=lambda u: H.SEQ.index(u) if u in H.SEQ else -1)
for _ in range(6):
    changed = False
    for u in order:
        i, kind = H.definer[u]
        old = v[u]
        E3._solvevar(v, ns, u, i, kind[0])
        if v[u] != old:
            changed = True
    if not changed:
        break
seedB = {f: v[f] for f in eng.FREE if v[f] != 0}
vB = eng.forward(seedB)
badB = eng.badatoms(vB)
FB = sorted(fscore.fails(badB))
print(f'BASELINE B (deliverable un-corrupted in place): score {fscore.score(badB)}, '
      f'{len(FB)} failures, {len(badB)} bad atoms')
print(f'  freed values put back: '
      f'{ {u: (len(str(abs(v[u]))) if v[u] else 0) for u in freed} } (digits)')
print(f'  demoted atoms now zero? '
      f'{ {a: (a not in badB) for a in demote} }')

# ---------- compare ----------
sA, sB = set(FA), set(FB)
print(f'\n=== COMPARISON ===')
print(f'  |A| {len(sA)}   |B| {len(sB)}   |A n B| {len(sA & sB)}')
print(f'  in A not B: {sorted(sA - sB)}')
print(f'  in B not A: {sorted(sB - sA)}')
print(f'  IDENTICAL: {sA == sB}')
print(f'\n  A = {FA}')
print(f'  B = {FB}')
print(f'\n  INTERSECTION ({len(sA & sB)}): {sorted(sA & sB)}')

json.dump({'A': FA, 'B': FB, 'intersection': sorted(sA & sB),
           'identical': sA == sB,
           'A_only': sorted(sA - sB), 'B_only': sorted(sB - sA)},
          open('baseline_sets.json', 'w'), indent=1)
print('\nwrote baseline_sets.json')
