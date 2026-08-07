import sys, os, json
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

def build(lock_extra=()):
    v=[0]*L.NVARS
    v[542]=1; v[438]=1
    fw.forward(v)
    v[14853]=v[12186]; v[16742]=v[24908]; v[8386]=0; v[21868]=0
    fw.forward(v)
    c0=L.polys[688][()]; m=8863713
    G0=(-c0*pow(m,-1,P))%P
    v[30213]=G0; v[22820]=0; v[7497]=(c0+m*G0)//P
    c0b=L.polys[1618][()]
    v[22162]=-c0b; v[14393]=0; v[11436]=0
    fw.forward(v)
    for a,t in [(13438,13153),(13440,20386),(36040,22917),(36042,26867)]:
        x=fw.solve_lin(a,t,v)
        if x is not None: v[t]=x
        fw.forward(v)
    locked=set([542,438,14853,16742,8386,21868,30213,22820,7497,22162,14393,11436,
                13153,20386,22917,26867])|set(lock_extra)
    for it in range(25):
        b=fw.bad_checks(v)
        if not b: break
        prog=False
        for a in b:
            if fw.evalpoly(L.polys[a],v)==0: continue
            for t in free_cands(a, locked):
                x=fw.solve_lin(a,t,v)
                if x is not None and x!=v[t]:
                    old=v[t]; v[t]=x; fw.forward(v)
                    if fw.evalpoly(L.polys[a],v)==0:
                        prog=True; break
                    v[t]=old; fw.forward(v)
        if not prog: break
    fw.forward(v)
    return v
v=build()
b=fw.bad_checks(v); av=L.all_atom_values(v); f=L.failing_eqs(av)
print(f"bad_checks={len(b)} failing={len(f)} score={L.NEQ-len(f)}")
print("bad:", b)
def fmt(a, lim=300):
    parts=[]
    for mm,c in sorted(L.polys[a].items(), key=lambda kv:(len(kv[0]),kv[0])):
        s=('%+d'%c) if (c not in (1,-1) or not mm) else ('+' if c==1 else '-')
        if mm: s+='*'.join('x%d'%u for u in mm)
        parts.append(s)
    return ' '.join(parts)[:lim]
for a in b:
    print(f"  a{a} eqs={len(L.atom2eq.get(a,{}))}: {fmt(a)}")
    print("     free:", [(u,NAT[u]) for u in sorted(L.avars[a]) if L.definer.get(u) is None])
json.dump({str(i):v[i] for i in range(L.NVARS)}, open('build4.json','w'))
