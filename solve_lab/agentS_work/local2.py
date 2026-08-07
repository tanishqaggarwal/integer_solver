import sys, json, collections
sys.path.insert(0,'.')
import common as C
import harness as H, engine as E
P=C.P
d=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
vD=[0]*E.NV
for k,val in d.items(): vD[int(k.split('_')[1])]=int(val)
ns={'v':vD,'__builtins__':{}}
BAD=[23616,23617,36659,36660,36661,36662,36663,36664]
VS=sorted(set().union(*[set(H.avars[a]) for a in BAD]))
print("local vars:",VS)
AT=sorted(set().union(*[set(H.occ[u]) for u in VS]))
print("atoms touching them:",len(AT))
for a in AT:
    r=eval(H.acodes[a],ns)
    print(f"  a{a} {'BAD' if r else '   '} : {H.atoms[a][:150]}")
print("\ncurrent values:")
for u in VS:
    x=vD[u]
    print(f"  x_{u} occ={len(H.occ[u])} bits={x.bit_length()} val={'0' if x==0 else (str(x)[:30]+'..' if abs(x)>10**30 else x)}")
