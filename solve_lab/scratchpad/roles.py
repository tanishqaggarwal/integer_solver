import json, sys
sys.path.insert(0, '/home/user/integer_solver/solve_lab/scratchpad')
import atomlib as A
p = A.p
BASE = '/home/user/integer_solver/solve_lab'

# free inputs from harness definition: gate outputs = targets in gates.jsonl
gate_out = set()
gate_def = {}   # target -> (rhs, vids)
with open(BASE + '/atoms/gates.jsonl') as f:
    for line in f:
        d = json.loads(line)
        gate_out.add(d['t'])
        gate_def.setdefault(d['t'], []).append((d['rhs'], tuple(d['vids'])))
freeinp = set(v for v in range(A.NVARS) if v not in gate_out)

v = A.load_json(BASE + '/best_agentA_39022.json')

vars_of_interest = [642, 2099, 7068, 4432, 19964, 28730, 28599, 17325, 17499, 9413,
                    23754, 26874, 6947]
print("VAR ROLES & VALUES at 39022:")
for x in vars_of_interest:
    role = 'FREE' if x in freeinp else 'GATE'
    ndef = len(gate_def.get(x, []))
    natoms = len(A.VAR_ATOMS[x])
    val = v[x]
    vp = val % p
    # is value == p or small?
    tag = ''
    if vp == p - 0: tag = ''
    if val == p: tag = ' [==p]'
    elif val == -p: tag = ' [==-p]'
    elif abs(val) < 10**12: tag = f' [small={val}]'
    defstr = ''
    if x in gate_def:
        defstr = ' ; '.join(r for r, _ in gate_def[x][:2])
    print(f"  x_{x}: {role}  in {natoms} atoms, {ndef} gatedefs  val%p={vp}{tag}")
    if defstr:
        print(f"        def: {defstr}")
