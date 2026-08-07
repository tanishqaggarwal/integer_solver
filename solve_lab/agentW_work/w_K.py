"""TASK 2: reproduce frame B and settle O's |K| = 34."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
import frameB
fr = frameB.Frame([642,28730,29854,31864])
W = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
v=[0]*frameB.NV
for k,val in W.items(): v[int(k[2:]) if k.startswith('x_') else int(k)]=int(val)
fv0={u:v[u] for u in fr.free if v[u]!=0}
st0=frameB.State(fr,fv0)
assert st0.score()==39026
A=37887
NZ=sorted(st0.nz())
print('nonzero check atoms (region):', NZ)
U=set()
for a in NZ: U |= set(fr.SUPV.get(a,[]))
C=set(fr.SUPV.get(A,[]))
print('is a37887 a check atom in frame B?', A in fr.csup)
print('|U| free inputs reaching a region atom =', len(U))
print('|C| free carriers of a37887 (S) =', len(C))
print('|U & C| =', len(U&C))
print('|U | C| = |K| =', len(U|C))
print('U =', sorted(U))
print('C =', sorted(C))
print('U&C =', sorted(U&C))
print('K =', sorted(U|C))
json.dump({'U':sorted(U),'C':sorted(C),'K':sorted(U|C),'NZ':NZ}, open('w_K.json','w'), indent=1)
# per-atom support sizes
for a in NZ:
    print('  atom %d support %d : %s' % (a, len(fr.SUPV.get(a,[])), sorted(fr.SUPV.get(a,[]))))
