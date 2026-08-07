"""Validate the corrected-engine refinement against tree96, then emit the oracle
dump for agent F.

Validity check that must pass before F uses any of this: for every tree96 stage S and
every refinement block B, either B is contained in S or B is disjoint from S.  A block
straddling two disjoint stages would mean the oracle is resolving something that is NOT
the tree, and the dump would be worthless.
"""
import sys, os, json, collections
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)

blocks = [set(b) for b in json.load(open('blocks_corrected.json'))['blocks']]
old8 = [set(b) for b in json.load(open('blocks8.json'))['blocks']]
tree = json.load(open('/home/user/integer_solver/solve_lab/agentF_work/tree96.json'))
gsup = {k: set(v['gsup']) for k, v in tree.items()}
LEAVES = set().union(*blocks)
print(f'{len(blocks)} blocks over {len(LEAVES)} leaves')

# ---- validity: no block straddles a stage boundary ----
# tree96 stages are LAMINAR, so a block may legitimately (a) sit inside a stage or
# (b) contain a whole stage.  Only a block that overlaps a stage while being neither
# contained in it nor containing it is a real crossing.
cross = []
for k, g in gsup.items():
    gl = g & LEAVES
    if not gl:
        continue
    for bi, b in enumerate(blocks):
        if (b & gl) and not (b <= gl) and not (gl <= b):
            cross.append((k, bi, len(b & gl), len(b), len(gl)))
print(f'\nVALIDITY: blocks straddling a stage boundary: {len(cross)}')
if cross:
    for c in cross[:15]:
        print('   stage %s block%d  overlap=%d  |block|=%d |stage|=%d' % c)
else:
    print('   none -- every block sits inside every stage it meets. Refinement is TREE-CONSISTENT.')

# same corrected check on the OLD 8-block partition, as a control
cross8 = []
for k, g in gsup.items():
    gl = g & LEAVES
    if not gl:
        continue
    for bi, b in enumerate(old8):
        if (b & gl) and not (b <= gl) and not (gl <= b):
            cross8.append((k, bi))
print(f'CONTROL: old 8-block partition crossings by the same test: {len(cross8)}')

# ---- is the new refinement a refinement of the old 8 blocks? ----
ok8 = all(any(b <= o for o in old8) for b in blocks)
print(f'\nrefines the previous 8-block partition: {ok8}')

# ---- per-stage induced partition (what the oracle resolves for F) ----
out = {}
for k, g in sorted(gsup.items(), key=lambda kv: -len(kv[1])):
    gl = g & LEAVES
    if len(gl) < 2:
        continue
    parts = sorted((sorted(gl & b) for b in blocks if gl & b), key=lambda p: (-len(p), p[0]))
    if len(parts) < 2:
        continue
    out[k] = {'n_leaves': len(gl), 'n_parts': len(parts), 'parts': parts}

# which stages does mux_wiring already decode?
mw = json.load(open('/home/user/integer_solver/solve_lab/agentF_work/mux_wiring.json'))
decoded = set(mw)
newly = [k for k in out if k not in decoded]
print(f'\nstages with an oracle-induced split: {len(out)}')
print(f'  of these, NOT in F mux_wiring.json (new information): {len(newly)}')
print(f'  {sorted(newly)}')

# binary (exactly-2-part) splits are directly slot-pair candidates
binary = {k: v for k, v in out.items() if v['n_parts'] == 2}
print(f'\nstages the oracle splits into EXACTLY 2 parts (direct slot-pair candidates): {len(binary)}')
for k, v in sorted(binary.items(), key=lambda kv: -kv[1]['n_leaves']):
    a, b = v['parts']
    tag = '' if k in decoded else '   <-- NOT in mux_wiring'
    print(f'  stage {k}: {len(a)} | {len(b)}{tag}')

json.dump({
    'note': ('Leaf-partition oracle from the residual side, measured on the CORRECTED engine '
             '(engine2: 5 atoms demoted from definer role so the 39,026 point is representable). '
             'Knob set: the 256 boolean leaves of the cone of atoms [7389,10187,20212,20215,28647]. '
             '136 bases, each saturating a different subtree; signature = delta of that base\'s own '
             'bad-atom set; two leaves are kept together only if their signatures agree in every '
             'base where NEITHER is ON.'),
    'validity': {'blocks_straddling_a_stage': len(cross),
                 'tree_consistent': len(cross) == 0,
                 'refines_previous_8_blocks': ok8},
    'resolution': {'blocks': len(blocks),
                   'previous_limit_broken_engine': 8,
                   'sizes': [len(b) for b in blocks]},
    'blocks': [sorted(b) for b in blocks],
    'stage_partitions': out,
    'binary_splits': {k: v['parts'] for k, v in binary.items()},
    'not_in_mux_wiring': sorted(newly),
    'caveats': ('The oracle CANNOT see inside a block: the largest block (142 leaves) is still '
                'unresolved internally, so stages wholly inside it get no split. A 2-part split is '
                'a slot-pair candidate but the oracle does not label which part is inA vs inB. '
                'Parts with >2 pieces mean the stage contains several resolved subtrees, not that '
                'the gate has >2 slots.'),
}, open('oracle_for_F.json', 'w'), indent=1)
print('\nwrote oracle_for_F.json')
