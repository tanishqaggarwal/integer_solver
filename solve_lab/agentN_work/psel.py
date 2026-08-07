"""SELECTOR AXIS: is the mod-p rank gap of the region response invariant under the one input
held fixed in all 16 detach states — the selector setting?

The 256 selectors are pure free inputs of Frame(POOL).  The deliverable has exactly TWO on:
{2081, 24601}.  Every configuration priced so far inherits that setting, so `p` has entered the
frame through one fixed set of constants.  Here the setting is varied STRUCTURALLY, using the
laminar block hierarchy recovered by seltree.py (a selector's position in the OR tree), and the
region response is re-measured.

Phase 1 (`python3 psel.py size`) is a SIZE PROBE only: state, score, |R|, knob count k.
Phase 2 (`python3 psel.py price <tag>...`) runs the exact pricing of pgap.py on the ones that fit.
"""
import os, sys, json, time
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(2000000)
import ev, optN
from optN import fr, FREE, FR0, BASEFV, atom_eqs, _bits
from frameB import State

SEL = sorted(int(k) for k in json.load(open(os.path.join(HERE, 'leafpins.json'))))
SELSET = set(SEL)
TREE = json.load(open(os.path.join(HERE, 'runs', 'seltree.json')))
BLOCKS = TREE['blocks']
WIT_ON = [s for s in SEL if BASEFV.get(s, 0)]


def min_block(members):
    """smallest hierarchy block containing all of `members` — the structural 'subtree' they
    share.  size 253 == they only meet at the root; None == they never meet."""
    ms = set(members)
    best = None
    for b in BLOCKS:
        if ms <= set(b['members']):
            if best is None or b['size'] < best['size']:
                best = b
    return best


def state_for(on):
    """State in Frame(POOL) with the witness free values and exactly `on` selectors set to 1."""
    fv = dict(BASEFV)
    ons = set(on)
    for s in SEL:
        if s in ons:
            fv[s] = 1
        else:
            fv.pop(s, None)
    return State(fr, fv)


def region_of(st):
    NZ = set(st.nz())
    R = set()
    for q in NZ:
        R |= atom_eqs[q]
    return NZ, R


def knobs_of(R):
    atoms_R = set()
    for e in R:
        for c, a in ev.eq_terms[e][2]:
            atoms_R.add(a)
    cands = set()
    for q in atoms_R:
        if q in fr.csup:
            cands.update(FR0[bb] for bb in _bits(fr.csup[q]))
    return sorted(y for y in cands if y in FREE)


# ---------------------------------------------------------------- configuration generator
def configs():
    """Structural variations of the live set, NOT numerical perturbations."""
    C = []
    C.append(('baseline_2081_24601', WIT_ON))
    # --- axis 1: keep cardinality 2, move where the pair sits in the tree ---------------
    # (a) both live leaves DEEP IN ONE SUBTREE (minimal common block as small as possible)
    small = [b for b in BLOCKS if b['size'] == 2]
    for i, b in enumerate(small[:3]):
        C.append(('pair_same_block2_%d' % i, b['members']))
    six = [b for b in BLOCKS if b['size'] == 6]
    if six:
        m = six[0]['members']
        C.append(('pair_same_block6', [m[0], m[1]]))
    # (b) both on ONE SIDE of the root but in different sub-subtrees
    big = sorted([b for b in BLOCKS if 60 <= b['size'] <= 100], key=lambda b: b['size'])
    if big:
        m = big[0]['members']
        C.append(('pair_same_side_far', [m[0], m[-1]]))
    # (c) OPPOSITE SIDES of the root
    d1 = sorted([b for b in BLOCKS if b['depth'] == 1 and b['size'] > 100],
                key=lambda b: b['size'])
    if len(d1) >= 2:
        a, b = set(d1[0]['members']), set(d1[-1]['members'])
        ao = sorted(a - b)
        bo = sorted(b - a)
        if ao and bo:
            C.append(('pair_opposite_root', [ao[0], bo[0]]))
    # (d) the live set drawn from a DIFFERENT SUBTREE than the witness's
    #     24601's own chain vs a disjoint chain
    own = min_block([24601])
    chain = sorted([b for b in BLOCKS if 24601 in b['members']], key=lambda b: b['size'])
    near = [b for b in chain if b['size'] in (6, 10)]
    if near:
        mm = [s for s in near[0]['members'] if s != 24601]
        C.append(('pair_in_24601_subtree', mm[:2]))
    excl = [b for b in BLOCKS if b['size'] in (6, 10) and 24601 not in b['members']
            and 2081 not in b['members']]
    if excl:
        C.append(('pair_foreign_subtree', excl[0]['members'][:2]))
        if len(excl) > 1:
            C.append(('pair_foreign_subtree2', excl[-1]['members'][:2]))
    # (e) 2081 is structurally exceptional (never joins the 253-block) — vary it alone
    C.append(('only_2081', [2081]))
    C.append(('only_24601', [24601]))
    C.append(('none_live', []))
    # --- axis 2: MORE LEAVES LIVE -------------------------------------------------------
    b10 = [b for b in BLOCKS if b['size'] == 10]
    b18 = [b for b in BLOCKS if b['size'] == 18]
    if b10:
        C.append(('block10_all_live', b10[0]['members']))
        C.append(('block10_half_live', b10[0]['members'][:5]))
    if b18:
        C.append(('block18_all_live', b18[0]['members']))
    C.append(('wit_plus_one', WIT_ON + [s for s in SEL if s not in WIT_ON][:1]))
    C.append(('wit_plus_two', WIT_ON + [s for s in SEL if s not in WIT_ON][:2]))
    C.append(('spread4', [SEL[0], SEL[64], SEL[128], SEL[192]]))
    C.append(('spread8', SEL[::32]))
    C.append(('first16', SEL[:16]))
    C.append(('all_live', SEL))
    return C


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'size'
    out = []
    t0 = time.time()
    for tag, on in configs():
        t1 = time.time()
        st = state_for(on)
        NZ, R = region_of(st)
        kb = knobs_of(R) if len(R) <= 4000 else None
        mb = min_block(on) if 1 < len(on) <= 20 else None
        rec = dict(tag=tag, live=sorted(on), nlive=len(on), score=st.score(),
                   fails=len(st.fails), nz=len(NZ), R=len(R),
                   knobs=(len(kb) if kb is not None else -1),
                   min_block=(mb['size'] if mb else None),
                   secs=round(time.time() - t1, 1))
        out.append(rec)
        print('%-24s live=%-4d score=%-6d |nz|=%-5d |R|=%-5d knobs=%-5d minblk=%-5s %.1fs'
              % (tag, rec['nlive'], rec['score'], rec['nz'], rec['R'], rec['knobs'],
                 rec['min_block'], rec['secs']), flush=True)
    json.dump(out, open(os.path.join(HERE, 'runs', 'psel_size.json'), 'w'), indent=1)
    print('total %.1fs -> runs/psel_size.json' % (time.time() - t0))


if __name__ == '__main__':
    main()
