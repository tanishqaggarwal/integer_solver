"""List the growth candidates that actually free a new private variable."""
import sys, json
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentO_work')
import regiongrow as G, harness as H

OD = '/home/user/integer_solver/solve_lab/agentO_work'
R0 = G.R0
P0 = set(G.private_vars(R0))
E0 = sorted({e for x in R0 for e in G.EQCO[x]})
cand = set()
for a in R0:
    for u in H.avars[a]:
        cand |= set(H.occ[u])
cand -= set(R0)
out = []
for a in sorted(cand):
    P = set(G.private_vars(R0 + [a]))
    Eq = sorted({e for x in R0 + [a] for e in G.EQCO[x]})
    dP = sorted(P - P0)
    if dP:
        out.append({'atom': a, 'frees': dP, 'nE': len(Eq), 'dE': len(Eq) - len(E0),
                    'src': H.atoms[a][:70]})
out.sort(key=lambda r: r['dE'])
json.dump(out, open(OD + '/growcand.json', 'w'), indent=1)
print('adjacent', len(cand), 'useful', len(out))
for r in out:
    print(f"  a{r['atom']}: frees {r['frees']} +{r['dE']} eqs (|E|={r['nE']})  {r['src']}")
