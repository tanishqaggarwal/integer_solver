import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
P=L.P
v=[int(x) for _,x in sorted(((int(k),val) for k,val in json.load(open('build2.json')).items()))]
fw.forward(v)
NAT={u:len(L.var_atoms[u]) for u in range(L.NVARS)}
def free_cands(a):
    out=[]
    for u in L.avars[a]:
        if L.definer.get(u) is not None: continue
        if any(mm.count(u)>1 for mm in L.polys[a]): continue
        out.append(u)
    out.sort(key=lambda u:(NAT[u],u))
    return out
# stage 1: the four simple load pins (their free data input)
for a,t in [(13438,13153),(13440,20386),(36040,22917),(36042,26867)]:
    x=fw.solve_lin(a,t,v)
    print(f"a{a}: solve x{t} -> {'OK' if x is not None else 'FAIL'} {str(x)[:30]}")
    if x is not None: v[t]=x
    fw.forward(v)
b=fw.bad_checks(v); print("after pins:", b)
# stage 2: big single-equation checks -> private slack var (fewest atoms)
for _ in range(8):
    b=fw.bad_checks(v)
    if not b: break
    prog=False
    for a in b:
        if fw.evalpoly(L.polys[a],v)==0: continue
        for t in free_cands(a):
            x=fw.solve_lin(a,t,v)
            if x is not None and x!=v[t]:
                old=v[t]; v[t]=x; fw.forward(v)
                if fw.evalpoly(L.polys[a],v)==0:
                    print(f"a{a}: fixed via x{t} (natoms={NAT[t]})"); prog=True; break
                v[t]=old; fw.forward(v)
    if not prog: break
fw.forward(v)
b=fw.bad_checks(v); av=L.all_atom_values(v); f=L.failing_eqs(av)
print(f"FINAL bad_checks={len(b)} failing={len(f)} score={L.NEQ-len(f)}")
print("bad:", b)
json.dump({str(i):v[i] for i in range(L.NVARS)}, open('build3.json','w'))
