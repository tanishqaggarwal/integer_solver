import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep, tri7
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); M=8640431
v=[0]*L.NVARS
for k,x in json.load(open(os.path.join(HERE,'data','quad3_hit.json'))).items(): v[int(k)]=int(x)
fw.forward(v)
LOCK={31339,33708,490,91,19750,7497,22820,14853,14393,11436,14515,16742,22162,30213,8386,21868,
      16441,28955,2751,18751}
print("start bad:", fw.bad_checks(v), "gamma%M:", (v[12000]//P)%M)
best=(len(L.failing_eqs(L.all_atom_values(v))), [x for x in v])
for rnd in range(8):
    tri7.close_all(v, LOCK, rounds=8, verbose=True)
    bad=fw.bad_checks(v); f=L.failing_eqs(L.all_atom_values(v))
    print(f" outer{rnd}: bad={len(bad)} failing={len(f)} score={L.NEQ-len(f)} gamma%M={(v[12000]//P)%M if v[12000]%P==0 else 'n/a'} {bad[:12]}", flush=True)
    if len(f)<best[0]: best=(len(f),[x for x in v])
    if not bad: break
print(f"BEST failing={best[0]} score={L.NEQ-best[0]}")
json.dump({str(i):best[1][i] for i in range(L.NVARS)}, open(os.path.join(HERE,'data','closehit2.json'),'w'))
