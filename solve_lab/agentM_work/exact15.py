"""Exact incident handle set, without the cofactor marker.

L's shortcut (equation e contains atom a <=> u_a in vars(e)) needs each residual atom to
have a unique free singly-occurring cofactor.  That holds for product-defined handles but
NOT for bare-defined ones (atom 36663 is literally `x_31864`) nor linearly-defined ones --
which is why my first pass returned 10 and dropped x31864, one of the deliverable's own
four.  I do not need the marker: H.eqt lists each equation's atoms directly, so

    handle h is incident  <=>  h's DEFINER ATOM appears in a target equation

is exact for every definer form.  rt(h) = how many target equations contain that atom.
"""
import sys, os, re, json, collections
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H
import engine as EB

BS = json.load(open('baseline_sets.json'))
T25 = BS['A']
LSTYLE = json.load(open('lcrit.json'))['L_style_baseline_fails']
D4 = [642, 28730, 29854, 31864]
NAMED = [23754, 35619, 9629, 37413, 34113, 28355]

ATOM2VAR = {H.definer[u][0]: u for u in H.SEQ}


def incident(targets, label):
    per = collections.Counter()
    for e in targets:
        for c, a in H.eqt[e][2]:
            if a >= 0 and c:
                per[a] += 1
    atoms = sorted(per)
    hs = {}
    for a in atoms:
        u = ATOM2VAR.get(a)
        if u is not None:
            hs[u] = (a, per[a])
    print(f'\n=== {label}: {len(targets)} equations, {len(atoms)} atoms, '
          f'{len(hs)} of those atoms are DEFINER atoms -> {len(hs)} incident handles ===')
    for u, (a, rt) in sorted(hs.items(), key=lambda kv: -kv[1][1]):
        t = H.atoms[a]
        form = ('product' if re.fullmatch(r'x_%d - x_\d+ \* x_\d+' % u, t)
                else 'bare' if t == 'x_%d' % u else 'linear/other')
        mark = '  <== DELIVERABLE' if u in D4 else ('  <== L-NAMED' if u in NAMED else '')
        print(f'  x{u:<6d} atom {a:<6d} rt {rt:2d}  [{form:12s}] {t[:44]}{mark}')
    return hs


h25 = incident(T25, 'against MY 25-equation baseline (the union L filtered on)')
h13 = incident(LSTYLE, 'against the L-style 13-equation baseline')

print('\n--- cross-checks ---')
print(f'deliverable four in the 25-set: { {u: (u in h25) for u in D4} }')
print(f'L-named handles in the 25-set : { {u: (u in h25) for u in NAMED} }')
pool32 = set(json.load(open('incident_pool.json'))['incident_handles'])
print(f'25-set vs my earlier 32-pool: |25-set| {len(h25)}, in 32 {len(set(h25)&pool32)}, '
      f'NOT in 32 {sorted(set(h25)-pool32)}')
print(f'13-set is subset of 25-set: {set(h13) <= set(h25)}')

order = sorted(h25.items(), key=lambda kv: -kv[1][1])
json.dump({'targets_25': T25, 'targets_13': LSTYLE,
           'incident_25': {str(u): {'atom': a, 'rt': rt} for u, (a, rt) in h25.items()},
           'incident_13': {str(u): {'atom': a, 'rt': rt} for u, (a, rt) in h13.items()},
           'ranked': [u for u, _ in order]},
          open('exact_incident.json', 'w'), indent=1)
print(f'\nwrote exact_incident.json  ({len(h25)} incident handles vs 25-eq baseline)')
