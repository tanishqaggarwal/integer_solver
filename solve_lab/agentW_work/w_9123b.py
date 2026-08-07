"""W: produce a CHECKER-VERIFIED assignment whose failing set contains 9123.
If it verifies, 9123 is reachable by an assignment, and the open half of the reconciliation is
purely 'U's configuration knob set vs my frame-B knob set' -- not 'is it reachable at all'."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
import w_setup2 as S, frameB
keep = [e for e in S.SAT if e != 9123]
sol = S.solve(keep + [12231])
print('integer solve for  (SAT \\ {9123}) + buy 12231 :', 'FOUND' if sol is not None else 'NONE')
assert sol is not None
ch = {u: S.fv0.get(u, 0) + sol.get(u, 0) for u in S.KNOB if sol.get(u, 0)}
st = S.st0.clone().set_free(ch)
print('exact score through frameB.State: %d   failing %s' % (st.score(), sorted(st.fails)))
print('knobs changed: %d  %s' % (len(ch), sorted(ch)))
out = 'w_pay9123_%d.json' % st.score()
json.dump({('x_%d' % i): str(st.v[i]) for i in range(frameB.NV) if st.v[i] != 0}, open(out, 'w'))
print('wrote', out)
