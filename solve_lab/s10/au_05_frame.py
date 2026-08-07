import os, sys, collections, random, json
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977
SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
E = [2554, 6816, 8124, 9123, 9421, 12231, 12270, 12350, 14584, 18673, 22044, 29125]
v = L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
av = L.all_atom_values(v)
VARS = [1329, 22665, 29854, 9118, 7075, 10903, 28961, 31864, 8731, 9413, 17499, 28730,
        7068, 2099, 642, 17325, 28599, 4432, 19964]
for u in VARS:
    d = L.definer.get(u)
    print(f'x_{u:<6} defined_by={d}  free={d is None}  val={str(v[u])[:60]}{"..." if len(str(v[u]))>60 else ""}  ==p:{v[u]==P}  bits={v[u].bit_length()}')
