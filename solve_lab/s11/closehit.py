import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep, tri7
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
v=[0]*L.NVARS
for k,x in json.load(open(os.path.join(HERE,'data','quad3_hit.json'))).items(): v[int(k)]=int(x)
fw.forward(v)
print("start bad:", fw.bad_checks(v), "failing:", len(L.failing_eqs(L.all_atom_values(v))))
print("  gamma%M =", (v[12000]//P)%8640431, " mirror0:", v[3719]%P==0 and v[25118]%P==0)
best=(len(L.failing_eqs(L.all_atom_values(v))), [x for x in v])
for rnd in range(6):
    tri7.close_all(v, set(), rounds=6, verbose=True)
    bad=fw.bad_checks(v); f=L.failing_eqs(L.all_atom_values(v))
    print(f" outer{rnd}: bad={len(bad)} failing={len(f)} score={L.NEQ-len(f)} {bad[:12]}", flush=True)
    if len(f)<best[0]: best=(len(f),[x for x in v])
    if not bad: break
print(f"BEST failing={best[0]} score={L.NEQ-best[0]}")
json.dump({str(i):best[1][i] for i in range(L.NVARS)}, open(os.path.join(HERE,'data','closehit.json'),'w'))
