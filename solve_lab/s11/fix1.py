"""Direct repair of the neighbourhood, then a gate ripple.

The three congruences local1 showed we need --  p | x9118, p | x8731, p | x28730 --  all sit on
variables we can move directly (x9118, x8731 are FREE; x28730 is defined by the very gate that
demands it).  So: round them down to multiples of p, set the four handles to the values the
checks want, ripple the gates, and count.
"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw, atomval
P=L.P
HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
src=sys.argv[1] if len(sys.argv)>1 else os.path.join(LAB,'best','new_instance_partial_39026.json')
v=load_raw(src)
def score(v):
    AV=[atomval(a,v) for a in range(L.NA)]
    F=L.failing_eqs(AV)
    return len(F),F,[a for a in range(L.NA) if AV[a]!=0]
n,F,NZ=score(v); print(f"start: failing={n} score={L.NEQ-n} nonzero atoms={len(NZ)}")

MODE=sys.argv[2] if len(sys.argv)>2 else 'round'
if MODE=='round':
    v[9118] -= v[9118] % P
    v[8731] -= v[8731] % P
elif MODE=='zero':
    v[9118] = 0
    v[8731] = 0
assert v[9118]%P==0 and v[8731]%P==0
# handles the checks want
v[1329]  = 5113045*v[7075]*v[9118] // P      # a35758/a35759 :  x29854 = p*x1329 = 5113045*x9118
v[29854] = 5113045*v[7075]*v[9118]
v[31864] = -v[7075]*v[8731]                  # a35761
v[10903] = v[31864] // P                     # a35760 :  x31864 = p*x10903
assert v[31864] == P*v[10903]
v[9413]  = v[28730] // P                     # a22230 :  x28730 = p*x9413
v[28730] = P*v[9413]
v[17325] = 0; v[642] = 0                     # a35762 :  x642 = p*x17325
v[7068]  = v[2099] + 7376877*v[642]          # a22229
v[21574] = 0; v[1844] = 0
AV=[atomval(a,v) for a in range(L.NA)]
print(" after direct set: nonzero atoms =", [a for a in range(L.NA) if AV[a]!=0])
n,F,NZ=score(v); print(f"   failing={n} score={L.NEQ-n}")

seeds={u:v[u] for u in (9118,8731,1329,29854,31864,10903,9413,28730,17325,642,7068,21574,1844)}
t0=time.time()
ch,st=L.ripple(v, seeds)
print(f" ripple: changed {len(ch)} vars in {st} steps ({time.time()-t0:.0f}s)")
n,F,NZ=score(v); print(f"   failing={n} score={L.NEQ-n} nonzero atoms={len(NZ)}")
print("   nonzero atoms:", NZ[:40])
print("   failing eqs:", F[:40])
json.dump({('x_%d'%i):v[i] for i in range(L.NVARS)}, open(os.path.join(HERE,'data',f'fix1_{MODE}.json'),'w'))
