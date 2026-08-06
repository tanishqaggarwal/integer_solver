import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw, atomval
P=L.P
HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
src=sys.argv[1] if len(sys.argv)>1 else os.path.join(LAB,'best','new_instance_partial_39026.json')
MODE=sys.argv[2] if len(sys.argv)>2 else 'round'
v=load_raw(src)
def rep(v,tag):
    AV=[atomval(a,v) for a in range(L.NA)]
    F=L.failing_eqs(AV); NZ=[a for a in range(L.NA) if AV[a]!=0]
    print(f"  {tag}: failing={len(F)} score={L.NEQ-len(F)} nonzero atoms={len(NZ)} {NZ[:30]}")
    return F,NZ
rep(v,'start')
n9 = v[9118]-(v[9118]%P) if MODE=='round' else 0
n8 = v[8731]-(v[8731]%P) if MODE=='round' else 0
seeds={9118:n9, 8731:n8,
       1329: 5113045*v[7075]*n9//P, 29854: 5113045*v[7075]*n9,
       31864: -v[7075]*n8, 10903: -v[7075]*n8//P,
       9413: v[28730]//P, 28730: P*(v[28730]//P),
       17325: 0, 642: 0, 7068: v[2099], 21574: 0, 1844: 0}
t0=time.time(); ch,st=L.ripple(v, seeds)
print(f"  ripple: changed {len(ch)} vars, {st} steps ({time.time()-t0:.0f}s)")
F,NZ=rep(v,'after ripple')
for a in NZ:
    print(f"    a{a} out={L.atom_out.get(a)} eqs={len(L.atom2eq.get(a,{}))}")
json.dump({('x_%d'%i):v[i] for i in range(L.NVARS)}, open(os.path.join(HERE,'data',f'fix2_{MODE}.json'),'w'))
