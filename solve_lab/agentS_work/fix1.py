import sys, json, collections
sys.path.insert(0,'.')
import common as C
import harness as H, engine as E
P=C.P
d=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
vD=[0]*E.NV
for k,val in d.items(): vD[int(k.split('_')[1])]=int(val)
def score(v,label):
    ns={'v':v,'__builtins__':{}}
    bad={i:eval(H.acodes[i],ns) for i in range(len(H.atoms))}
    bad={i:r for i,r in bad.items() if r}
    ff=E.eqfails(bad)
    print(f"{label}: atoms={sorted(bad)} fails={len(ff)} SCORE={39033-len(ff)}  eqs={ff}")
    return bad,ff
score(vD,"deliverable")
print("\nmod-p status of the handle numerators:")
for u in (642,28730,29854,31864):
    print(f"  x_{u} % p == 0 ? {vD[u]%P==0}   (bits {vD[u].bit_length()})")
v=list(vD); v[31864]=0; v[10903]=0
score(v,"+ x_31864=0, x_10903=0")
