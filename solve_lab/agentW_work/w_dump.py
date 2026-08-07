"""Dump the 'buy eq12231, break eq2554' trade as a full assignment for checker.py."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w_setup as S
import frameB
sol = S.solve([e for e in S.SAT if e != 2554] + [12231])
assert sol is not None
ch = {u: S.fv0.get(u, 0) + sol.get(u, 0) for u in S.KNOB if sol.get(u, 0)}
st = S.st0.clone().set_free(ch)
print('model says score %d failing %s' % (st.score(), sorted(st.fails)))
json.dump({('x_%d' % i): str(st.v[i]) for i in range(frameB.NV) if st.v[i] != 0},
          open('w_trade_12231_break2554.json', 'w'))
print('knob values changed vs witness:', {u: (S.fv0.get(u,0), ch[u]) for u in list(ch)[:6]})
print('vars differing from witness:', sum(1 for i in range(frameB.NV) if st.v[i] != S.st0.v[i]))
