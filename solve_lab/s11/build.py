import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw
v=[0]*L.NVARS
v[542]=1; v[438]=1
fw.forward(v)
print("after bits: x12186=%s x24908=%s x15298=%d x15574=%d"%(v[12186],v[24908],v[15298],v[15574]))
# kill the core: n = x14853-x12186 = 0, m = x24908-x16742 = 0
v[14853]=v[12186]; v[16742]=v[24908]
v[8386]=0; v[21868]=0
fw.forward(v)
print("n=%s m=%s x35389=%s x6671=%s"%(v[29322],v[3558],v[35389],v[6671]))
b=fw.bad_checks(v); av=L.all_atom_values(v); f=L.failing_eqs(av)
print(f"bad_checks={len(b)} failing={len(f)} score={L.NEQ-len(f)}")
print("bad:", b)
def fmt(a, lim=260):
    parts=[]
    for m,c in sorted(L.polys[a].items(), key=lambda kv:(len(kv[0]),kv[0])):
        s=('%+d'%c) if (c not in (1,-1) or not m) else ('+' if c==1 else '-')
        if m: s+='*'.join('x%d'%u for u in m)
        parts.append(s)
    return ' '.join(parts)[:lim]
for a in b:
    print(f"  a{a} eqs={len(L.atom2eq.get(a,{}))}: {fmt(a)}")
    for u in sorted(L.avars[a]):
        if L.definer.get(u) is None:
            print(f"       FREE x{u} (in {len(L.var_atoms[u])} atoms)")
json.dump({str(i):v[i] for i in range(L.NVARS)}, open('build0.json','w'))
