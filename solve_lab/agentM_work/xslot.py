"""Cross-check the corrected-engine refinement against F's OWN decode (mux_wiring.json).

A refinement block that cuts a stage is only legitimate if the cut is the stage's own
slot boundary.  Slot leaf-support is computed from F's mux_wiring selector variables ->
defining atom's cone -> intersected with the 256 leaves, exactly as in the earlier
root / 19538 / 10649 confirmations.  No use of my own measurement in building them.
"""
import sys, os, json, re, collections
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine as E
import mcore2 as M

F_DIR = '/home/user/integer_solver/solve_lab/agentF_work'
mw = json.load(open(os.path.join(F_DIR, 'mux_wiring.json')))
tree = json.load(open(os.path.join(F_DIR, 'tree96.json')))
gsup = {k: set(v['gsup']) for k, v in tree.items()}
LEAVES = set(M.bools())
blocks = [set(b) for b in json.load(open('blocks_corrected.json'))['blocks']]
old8 = [set(b) for b in json.load(open('blocks8.json'))['blocks']]

VARRE = re.compile(r'x(\d+)')


def slot_support(entry):
    """Leaf support of one slot, from its selector variables' cones."""
    sup = set()
    for sel, gates in entry:
        ids = {int(sel)}
        for _, expr in gates:
            ids |= {int(m) for m in VARRE.findall(expr)}
        for u in ids:
            for a in H.occ.get(u, ()):
                try:
                    sup |= set(E.cone(a)[1])
                except Exception:
                    pass
    return sup & LEAVES


print('stage | slotA | slotB (from F mux_wiring) | my induced parts | verdict')
agree = dis = 0
report = {}
for k, ent in mw.items():
    if 'inA' not in ent or 'inB' not in ent:
        continue
    gl = gsup.get(k, set()) & LEAVES
    if len(gl) < 2:
        continue
    A = slot_support(ent['inA']) & gl
    B = slot_support(ent['inB']) & gl
    parts = sorted((gl & b for b in blocks if gl & b), key=lambda s: -len(s))
    if len(parts) < 2:
        continue
    # does my partition's coarsest 2-way grouping match {A,B}?
    match = None
    if A and B and not (A & B) and (A | B) == gl:
        # each of my parts must lie wholly in A or wholly in B
        ok = all((p <= A) or (p <= B) for p in parts)
        match = 'CONSISTENT with slot pair' if ok else 'CUTS ACROSS the slot pair'
    else:
        match = f'slot supports unusable (|A|={len(A)},|B|={len(B)},|A&B|={len(A&B)},cover={len(A|B)}/{len(gl)})'
    if match.startswith('CONSISTENT'):
        agree += 1
    elif match.startswith('CUTS'):
        dis += 1
    report[k] = {'nleaves': len(gl), 'slotA': len(A), 'slotB': len(B),
                 'my_parts': [len(p) for p in parts], 'verdict': match}
    print(f'  {k}: {len(A)}|{len(B)}  mine {[len(p) for p in parts]}  -> {match}')

print(f'\nCONSISTENT: {agree}   CUTS ACROSS: {dis}')

# same test for the old 8-block partition, as a control
agree8 = dis8 = 0
for k, ent in mw.items():
    if 'inA' not in ent or 'inB' not in ent:
        continue
    gl = gsup.get(k, set()) & LEAVES
    if len(gl) < 2:
        continue
    A = slot_support(ent['inA']) & gl
    B = slot_support(ent['inB']) & gl
    parts = [gl & b for b in old8 if gl & b]
    if len(parts) < 2 or not (A and B and not (A & B) and (A | B) == gl):
        continue
    if all((p <= A) or (p <= B) for p in parts):
        agree8 += 1
    else:
        dis8 += 1
print(f'CONTROL old-8: CONSISTENT {agree8}, CUTS ACROSS {dis8}')

json.dump(report, open('xslot_report.json', 'w'), indent=1)
print('wrote xslot_report.json')
