"""AUDIT of O's collateral accounting: 'every one of the 7 is individually buyable and every
purchase costs exactly eq8680'.  Solve each single purchase against EACH essential break row,
then price the result EXACTLY through frameB.State (not through the linear model)."""
import sys, os, itertools, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w_setup as S
import frameB

ESS = [2554, 6816, 8124, 9123, 9421, 'S']
out = {}
for f in S.FAIL:
    row = []
    for r in ESS:
        keep = [e for e in S.SAT if e != r]
        sol = S.solve(keep + [f])
        if sol is None:
            row.append((str(r), None)); continue
        ch = {u: S.fv0.get(u, 0) + sol.get(u, 0) for u in S.KNOB if sol.get(u, 0)}
        st = S.st0.clone().set_free(ch)
        row.append((str(r), st.score(), sorted(st.fails)))
    out[f] = row
    print('buy eq%d:' % f, flush=True)
    for e in row:
        if e[1] is None: print('    break %-6s : model-infeasible' % e[0])
        else: print('    break %-6s : EXACT score %d  failing %s' % (e[0], e[1], e[2]))
json.dump({str(k): [[str(x) for x in e] for e in v] for k, v in out.items()},
          open('w_audit.json', 'w'), indent=1)
