"""S12 step 12: genuine SECOND-ORDER synergy in the exhaustive pair data?
A pair is second-order only if knobs(u,z) > max(knobs(u), knobs(z))."""
import os, sys, json, collections
HERE='/home/user/integer_solver/solve_lab/s10'
D = json.load(open(os.path.join(HERE,'ac_single.json')))
S = {int(z): r[0]['knobs'] for z, r in D['res'].items()}
C = {int(z): r[0]['lost'] for z, r in D['res'].items()}
res = json.load(open(os.path.join(HERE,'ac_pairs.json')))['res']
syn = []
for k, lost, chk, atoms, u, z, sc in res:
    m = max(S[u], S[z])
    if k > m: syn.append((k - m, k, m, lost, u, z, S[u], S[z], C[u], C[z], sc))
syn.sort(key=lambda t: (-t[0], t[3]))
print(f'pairs measured: {len(res)} grew the support; genuinely second-order '
      f'(pair beats both singles): {len(syn)}')
h = collections.Counter(t[0] for t in syn)
print('synergy histogram (extra knobs over the best single):', dict(sorted(h.items())))
print('\ncheapest genuine second-order pairs '
      '(extra, pair_knobs, best_single, eqs_lost, x_u, x_z, k_u, k_z, cost_u, cost_z, score):')
for t in syn[:25]: print(f'  {t}')
both0 = [t for t in syn if t[6] == 0 and t[7] == 0]
print(f'\npairs where BOTH singles give zero knobs but the pair gives knobs: {len(both0)}')
for t in both0[:20]: print(f'  {t}')
