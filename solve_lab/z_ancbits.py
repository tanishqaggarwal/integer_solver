import os,sys,json
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
p=H.p
pins=json.load(open('pinrec.json'))
selectors=set(r[1] for r in pins)
FAILS11=[2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125]
RIPPLE16=[697,1985,5225,10815,16048,17784,17801,22402,23667,24721,27124,28737,29638,29959,35935,37431]
def bitanc(eqs,label):
    allvars=set()
    for i in eqs: allvars|=H.eqvars[i]
    anc=set()
    for v in allvars: anc|=H.anc.get(v,{v})
    bits=anc & selectors
    print(f"{label}: {len(eqs)} eqs, {len(allvars)} vars, {len(anc)} free ancestors, {len(bits)} selector-bits in ancestry")
    print(f"   bits: {sorted(bits)}")
    return bits
b11=bitanc(FAILS11,"11 fails")
b16=bitanc(RIPPLE16,"16 ripple")
b27=bitanc(FAILS11+RIPPLE16,"27 combined")
# Also: full closure including gap-knob descendants. The ripple only appears AFTER closing gaps (moving x_7068,x_4432).
# So include eqs downstream of x_7068,x_4432. Get their descendant gate outputs, then eqs using them.
from collections import defaultdict
desc=set()
for t in H.order:
    if H.anc[t]&{4432,7068}: desc.add(t)
desc|={4432,7068}
downstream_eqs=set()
for i,vs in enumerate(H.eqvars):
    if vs & desc: downstream_eqs.add(i)
print(f"\neqs downstream of x_4432/x_7068: {len(downstream_eqs)}")
allvars=set()
for i in downstream_eqs: allvars|=H.eqvars[i]
anc=set()
for v in allvars: anc|=H.anc.get(v,{v})
bits=anc&selectors
print(f"   their ancestry: {len(anc)} free, {len(bits)} bits: {sorted(bits)}")
