import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
v=[0]*L.NVARS
for b in (542,47,438,91): v[b]=1
fw.forward(v)
print("a=%d b=%d c=%d d=%d U=%d V=%d x15298=%d"%(v[8599],v[21839],v[7304],v[25956],v[7715],v[34554],v[15298]))
print("x38170(a*b)=%d  x3896(c*d)=%d"%(v[38170],v[3896]))
bad=fw.bad_checks(v); av=L.all_atom_values(v); f=L.failing_eqs(av)
print(f"bad_checks={len(bad)} failing={len(f)} score={L.NEQ-len(f)}")
print("bad:", bad)
# controllability test
for tgt,handle in [(12186,5096),(1308,14515),(24908,19750),(19083,21589)]:
    old=v[handle]; before=v[tgt]
    v[handle]=old+1; fw.forward(v)
    print(f"  x{tgt}: d/dx{handle} = {v[tgt]-before}")
    v[handle]=old; fw.forward(v)
