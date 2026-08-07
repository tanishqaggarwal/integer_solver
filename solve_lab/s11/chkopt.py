import sys, os, json, itertools, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep
from zsolve import solve_int
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(300000)
LAB=os.path.join(HERE,'..')
v=[0]*L.NVARS
for k,x in json.load(open(os.path.join(LAB,'best','new_instance_partial_39026.json'))).items():
    v[int(k[2:]) if k.startswith('x_') else int(k)]=int(x)
fw.forward(v)
av=L.all_atom_values(v)
S=sorted(L.failing_eqs(av))
print("checkpoint failing:", len(S), S)
region=set()
for e in S: region |= set(L.eq_atoms[e][2])
knobs=[a for a in region if set(L.atom2eq.get(a,{})).issubset(set(S))]
print("region atoms:", len(region), " knobs:", sorted(knobs))
steps=[]
for a in knobs:
    try: hs,_=deep.handles(v,a,locked=set())
    except Exception: hs=[]
    pr=[(t,d) for t,d in hs if len(L.var_atoms[t])==1 and d]
    tag = pr[0] if pr else None
    print(f"   a{a}: handles={[(t,len(L.var_atoms[t])) for t,_ in hs][:6]} private={[t for t,_ in pr]}")
    if pr: steps.append((a,pr[0][0],pr[0][1]))
print("realisable knobs:", [(a,t) for a,t,_ in steps])
if steps:
    rows=[];rhs=[]
    for e in S:
        mult,sq,co=L.eq_atoms[e]
        rows.append([co.get(a,0)*d for a,t,d in steps])
        rhs.append(-sum(c*av[a] for a,c in co.items()))
    for drop in range(0,len(S)+1):
        found=None
        for combo in itertools.combinations(range(len(S)),drop):
            keep=[i for i in range(len(S)) if i not in combo]
            x=solve_int([rows[i] for i in keep],[rhs[i] for i in keep])
            if x is not None: found=(combo,x); break
        if found:
            print(f"  MAX SATISFIABLE = {len(S)-drop} of {len(S)}  -> failing {drop}, score {L.NEQ-drop}")
            break
    else:
        print("  none satisfiable")
