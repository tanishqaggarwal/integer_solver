"""Which vars do the 8 nonzero atoms define in E's orientation?"""
import sys, os, json, collections
os.chdir('/home/user/integer_solver/solve_lab/agentM_work')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H, engine as E

BAD8 = [23616, 23617, 36659, 36660, 36661, 36662, 36663, 36664]
# invert definer: atom -> var
atom2var = {}
for u in H.SEQ:
    i, kind = H.definer[u]
    atom2var.setdefault(i, []).append(u)

print('atom -> defined var, for the 8 nonzero atoms of the deliverable:')
for a in BAD8:
    print(f'  atom {a}: defines {atom2var.get(a)}   (is a definer: {a in atom2var})')

print()
print('total distinct definer atoms:', len(atom2var))
print('atoms used to define >1 var:', sum(1 for a, l in atom2var.items() if len(l) > 1))

# the deliverable's vector
def load_vec(path):
    d = json.load(open(path))
    v = [0] * H.NV
    for k, val in d.items():
        v[int(k.split('_')[1])] = int(val)
    return v
vd = load_vec('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')

ns = {'v': vd, '__builtins__': {}}
print()
print('value of each of the 8 atoms AT THE DELIVERABLE:')
for a in BAD8:
    print(f'  atom {a}: {str(eval(H.acodes[a], ns))[:40]}')

# and at E's forward
import pickle
d1 = pickle.load(open('diag1.pkl', 'rb'))
vf = E.forward(d1['seed'])
ns2 = {'v': vf, '__builtins__': {}}
print()
print("value of each of the 8 atoms AT E's FORWARD:")
for a in BAD8:
    print(f'  atom {a}: {str(eval(H.acodes[a], ns2))[:40]}')

print()
print('atom text for the 8:')
for a in BAD8:
    print(f'  {a}: {H.atoms[a][:160]}')
