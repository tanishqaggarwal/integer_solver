import sys, json, collections
sys.path.insert(0,'.')
import common as C
import harness as H, engine as E
P=C.P
d=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
vD=[0]*E.NV
for k,val in d.items(): vD[int(k.split('_')[1])]=int(val)
ns={'v':vD,'__builtins__':{}}
def val(u): return vD[u]
for u in [17499,22665,28599,28961]:
    print(f"x_{u} == p ? {vD[u]==P}")
KEY=[642,2099,7068,7075,8731,9118,28730,29854,31864,10903,1329,17325,9413,21279,25297,37158]
for u in KEY:
    print(f"\nx_{u} occ={len(H.occ[u])} val_bits={vD[u].bit_length()}:")
    for a in H.occ[u]:
        r=eval(H.acodes[a],ns)
        print(f"    a{a} {'BAD' if r else 'ok '}: {H.atoms[a][:130]}")
