import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, close2
P=L.P
BITS=(542,47,438,91)
C0=L.polys[688][()]; MM=8863713; G0=(-C0*pow(MM,-1,P))%P; C0B=L.polys[1618][()]
TH={int(k):x for k,x in json.load(open('theta_solveB.json')).items()}
v=[0]*L.NVARS
for b in BITS: v[b]=1
for k,x in TH.items(): v[k]=x
fw.forward(v)
v[14853]=v[12186]; v[16742]=v[24908]
v[30213]=G0; v[22820]=0; v[7497]=(C0+MM*G0)//P
v[22162]=-C0B; v[14393]=0; v[11436]=0
fw.forward(v)
LOCK=set(BITS)|set(TH)|{14853,16742,30213,22820,7497,22162,14393,11436}
v,b=close2.close(v, LOCK)
av=L.all_atom_values(v); f=L.failing_eqs(av)
print(f"FINAL bad={len(b)} failing={len(f)} score={L.NEQ-len(f)}")
print("bad:", b)
json.dump({str(i):v[i] for i in range(L.NVARS)}, open('run2.json','w'))
