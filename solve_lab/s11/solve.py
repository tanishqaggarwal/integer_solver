"""Full constructive solve."""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
P=L.P
NAT={u:len(L.var_atoms[u]) for u in range(L.NVARS)}
def free_cands(a, locked):
    out=[]
    for u in L.avars[a]:
        if L.definer.get(u) is not None or u in locked: continue
        if any(mm.count(u)>1 for mm in L.polys[a]): continue
        out.append(u)
    out.sort(key=lambda u:(NAT[u],u)); return out

v=[0]*L.NVARS
BITS=(542,47,438,91)
for b in BITS: v[b]=1
fw.forward(v)

def structural():
    # core: n = x14853-x12186 = 0 ; m = x24908-x16742 = 0
    v[14853]=v[12186]; v[16742]=v[24908]
    fw.forward(v)
    # a29539 : x1308 == x14853   (via handle x14515, d/d = 1)
    v[14515]+= v[14853]-v[1308]
    # a26731 : x19083 == x16742  (via handle x21589, d/d = 1)
    v[21589]+= v[16742]-v[19083]
    v[8386]=0; v[21868]=0
    fw.forward(v)
    # arithmetic cluster
    c0=L.polys[688][()]; mm=8863713
    G0=(-c0*pow(mm,-1,P))%P
    v[30213]=G0; v[22820]=0; v[7497]=(c0+mm*G0)//P
    v[22162]=-L.polys[1618][()]; v[14393]=0; v[11436]=0
    fw.forward(v)

structural()
print("n=%d m=%d  x1308-x14853=%d  x19083-x16742=%d"%(v[29322],v[3558],v[1308]-v[14853],v[19083]-v[16742]))
b=fw.bad_checks(v); print("after structural:", len(b), b)

LOCK=set(BITS)|{14853,16742,30213,22162,5096,19750,14515,21589,7497,22820,14393,11436,8386,21868,12186,24908}
t0=time.time()
for outer in range(6):
    for it in range(40):
        b=fw.bad_checks(v)
        if not b: break
        prog=False
        for a in b:
            if fw.evalpoly(L.polys[a],v)==0: continue
            for t in free_cands(a, LOCK):
                x=fw.solve_lin(a,t,v)
                if x is not None and x!=v[t]:
                    old=v[t]; v[t]=x; fw.forward(v)
                    if fw.evalpoly(L.polys[a],v)==0: prog=True; break
                    v[t]=old; fw.forward(v)
        if not prog: break
    b=fw.bad_checks(v)
    av=L.all_atom_values(v); f=L.failing_eqs(av)
    print(f"outer {outer}: bad={len(b)} failing={len(f)} score={L.NEQ-len(f)} ({time.time()-t0:.0f}s) {b[:12]}")
    if not b: break
    structural()
fw.forward(v)
b=fw.bad_checks(v); av=L.all_atom_values(v); f=L.failing_eqs(av)
print(f"FINAL bad={len(b)} failing={len(f)} score={L.NEQ-len(f)}")
print("bad:", b)
json.dump({str(i):v[i] for i in range(L.NVARS)}, open('solve.json','w'))
