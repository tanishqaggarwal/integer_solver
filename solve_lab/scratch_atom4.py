import json
p=2**256-2**32-977
atoms=[]; reprs=[]; eqs=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f:
        dd=json.loads(line)
        atoms.append([(tuple(m),c) for m,c in dd['poly']])
        reprs.append(dd.get('repr',''))
        eqs.append(dd.get('eqs',[]))
for i in [18081,18084,29377,35321]:
    print(f"--- atom {i} ---")
    print("repr:", reprs[i])
    print("in eqs:", eqs[i])
    # variables
    vs=set()
    for m,c in atoms[i]:
        vs.update(m)
    print("vars:", sorted(vs))
