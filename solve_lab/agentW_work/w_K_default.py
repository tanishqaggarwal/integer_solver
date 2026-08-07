"""Same K construction in the DEFAULT orientation (no detach) -- reproduce T's 12/11/23/0."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
import frameB
for DET in ([], [642,28730,29854,31864]):
    fr = frameB.Frame(DET)
    W = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
    v=[0]*frameB.NV
    for k,val in W.items(): v[int(k[2:]) if k.startswith('x_') else int(k)]=int(val)
    fv0={u:v[u] for u in fr.free if v[u]!=0}
    st0=frameB.State(fr,fv0)
    A=37887
    NZ=sorted(st0.nz())
    U=set()
    for a in NZ: U |= set(fr.SUPV.get(a,[]))
    C=set(fr.SUPV.get(A,[]))
    print('DET=%s  free=%d  score=%d  nz=%s' % (DET, len(fr.free), st0.score(), NZ))
    print('   |U|=%d |C|=%d |U&C|=%d |K|=%d' % (len(U),len(C),len(U&C),len(U|C)))
    print('   U=%s' % sorted(U))
    print('   C=%s' % sorted(C))
