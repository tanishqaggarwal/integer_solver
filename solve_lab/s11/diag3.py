import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
v=[0]*L.NVARS
for k,x in json.load(open(os.path.join(HERE,'data','finish3.json'))).items(): v[int(k)]=int(x)
fw.forward(v)
print("bad:", fw.bad_checks(v), "failing:", len(L.failing_eqs(L.all_atom_values(v))))
print("x3719 %P:", v[3719]%P != 0, " x25118 %P:", v[25118]%P != 0, " x3896:", v[3896])
for a in [26719,26721,26723]:
    val=fw.evalpoly(L.polys[a],v)
    hs,base=deep.handles(v,a,locked=set())
    print(f"a{a}: val bits={val.bit_length()} val%P==0:{val%P==0}")
    for t,d in hs[:4]:
        print(f"    handle x{t} delta/P = {d//P if d%P==0 else 'not P-mult'} ; base%delta==0: {base%d==0}")
    print(f"    vars: {sorted(L.avars[a])}")
print("x12000 =", str(v[12000])[:40], " %P==0:", v[12000]%P==0)
print("x12926 %P==0:", v[12926]%P==0, " x21364 %P==0:", v[21364]%P==0)
