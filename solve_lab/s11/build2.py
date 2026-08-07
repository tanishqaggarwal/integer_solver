import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
P=L.P
v=[0]*L.NVARS
v[542]=1; v[438]=1
fw.forward(v)
v[14853]=v[12186]; v[16742]=v[24908]; v[8386]=0; v[21868]=0
fw.forward(v)
assert v[15298]==1 and v[5647]==0 and v[34606]==0, (v[15298],v[5647],v[34606])

# ---- a688 / a40608 : x14257 - 8863713*x18956 = c0 ----
c0 = L.polys[688][()]                     # a688 = c0 - x14257 + 8863713 x18956
m  = 8863713
G0 = (-c0 * pow(m, -1, P)) % P            # required x37892 mod p
num = c0 + m*G0
assert num % P == 0
Q = num // P
v[30213] = G0            # x37892 = x30213
v[22820] = 0
v[7497]  = Q             # x14257 = P*x7497
# ---- a1618 : x24468 - x32989 = -c0b ----
c0b = L.polys[1618][()]
v[22162] = -c0b          # x13682 = x22162
v[14393] = 0; v[11436] = 0
fw.forward(v)
print("x37892=%s"%str(v[37892])[:40], " x13682 ok:", v[13682]==-c0b)
print("a688 =", fw.evalpoly(L.polys[688],v))
print("a1618=", fw.evalpoly(L.polys[1618],v))
print("a40608=", fw.evalpoly(L.polys[40608],v))
b=fw.bad_checks(v); av=L.all_atom_values(v); f=L.failing_eqs(av)
print(f"bad_checks={len(b)} failing={len(f)} score={L.NEQ-len(f)}")
print("bad:", b)
json.dump({str(i):v[i] for i in range(L.NVARS)}, open('build2.json','w'))
