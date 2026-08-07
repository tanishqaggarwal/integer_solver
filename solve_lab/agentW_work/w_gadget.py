"""Pull one real stage gadget verbatim out of the instance and confirm its two congruences."""
import sys, os, json, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
import model
d = model.get(); atom_src = d['atom_src']; atom_vars = d['atom_vars']; eq_terms = d['eq_terms']
st = json.load(open('/home/user/integer_solver/solve_lab/agentQ_work/qstages.json'))['stages']
g = st[0]
print('stage 0 wires:', g)
want = {v for k, v in g.items() if k != 'kind'}
byvar = {}
for a, vs in enumerate(atom_vars):
    for v in vs: byvar.setdefault(v, []).append(a)
seen = set()
for name, v in sorted(g.items(), key=lambda t: str(t[0])):
    if name == 'kind': continue
    for a in byvar.get(v, []):
        if a in seen: continue
        s = atom_src[a]
        if len(s) < 260 and set(int(m) for m in re.findall(r'x_(\d+)', s)) & want:
            seen.add(a)
            print('  a%-6d %s' % (a, s))
