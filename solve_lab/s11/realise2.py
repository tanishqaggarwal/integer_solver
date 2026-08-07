import sys, os, json, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep
from zsolve import solve_int
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
v=[0]*L.NVARS
for k,x in json.load(open(os.path.join(HERE,'data','finish3.json'))).items(): v[int(k)]=int(x)
fw.forward(v)
av=L.all_atom_values(v)
S=sorted(L.failing_eqs(av))
def inner(vv):
    a2=L.all_atom_values(vv)
    return [sum(c*a2[a] for a,c in L.eq_atoms[e][2].items()) for e in S]
base=inner(v)
print("failing eqs:", len(S), " inner sums nonzero:", sum(1 for x in base if x), flush=True)
cands=set()
for a in [26719,26720,26721,26722,26723,28437]:
    try:
        hs,_=deep.handles(v,a,locked=set())
        for t,d in hs: cands.add(t)
    except Exception: pass
region=set()
for e in S: region |= set(L.eq_atoms[e][2])
for a in region:
    for u in L.avars[a]:
        if L.definer.get(u) is None: cands.add(u)
cands=sorted(cands)
print("candidate handles:", len(cands), flush=True)
cols=[];used=[]
for t in cands:
    old=v[t]
    v[t]=old+1; fw.forward(v); e1=inner(v)
    v[t]=old+2; fw.forward(v); e2=inner(v)
    v[t]=old; fw.forward(v)
    d1=[e1[i]-base[i] for i in range(len(S))]
    d2=[e2[i]-base[i] for i in range(len(S))]
    if all(d2[i]==2*d1[i] for i in range(len(S))) and any(d1):
        cols.append(d1); used.append(t)
print("EXACT LINEAR handles:", len(used), used, flush=True)
if used:
    M=[[cols[j][i] for j in range(len(used))] for i in range(len(S))]
    rhs=[-base[i] for i in range(len(S))]
    x=solve_int(M,rhs)
    print("integer solution:", "FOUND" if x else "NONE", flush=True)
    if x:
        for j,t in enumerate(used): v[t]+=x[j]
        fw.forward(v)
        f=L.failing_eqs(L.all_atom_values(v)); b=fw.bad_checks(v)
        print(f"AFTER: failing={len(f)} score={L.NEQ-len(f)} bad_checks={len(b)} {b[:10]}")
        json.dump({('x_%d'%i):v[i] for i in range(L.NVARS)}, open(os.path.join(HERE,'data','realise_named.json'),'w'))
    else:
        for drop in range(1,12):
            hit=None
            for combo in itertools.combinations(range(len(S)),drop):
                keep=[i for i in range(len(S)) if i not in combo]
                if solve_int([M[i] for i in keep],[rhs[i] for i in keep]) is not None:
                    hit=combo; break
            if hit:
                print(f"  recover {len(S)-drop} of {len(S)} -> failing {drop}, score {L.NEQ-drop}")
                break
