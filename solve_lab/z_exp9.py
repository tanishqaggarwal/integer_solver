import os,sys,json
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
atoms=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f: atoms.append(json.loads(line))
for ai in [42669,44342,45677,7450,7452]:
    a=atoms[ai]
    print(f"=== atom {ai}  (n_eq={a.get('n_eq')}, eqs={sorted(set(a['eqs']))[:3]}) ===")
    print(a['repr'])
    print()
