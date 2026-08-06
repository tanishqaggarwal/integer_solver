import sys, os, json, itertools, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, deep
from zsolve import solve_int
sys.set_int_max_str_digits(200000)
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
v=[0]*L.NVARS
for k,x in json.load(open(os.path.join(HERE,'data','finish3.json'))).items(): v[int(k)]=int(x)
fw.forward(v)
av=L.all_atom_values(v)
S=sorted(L.failing_eqs(av))
# candidate handles: exact-linear effect on S
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
# E = every equation any candidate can touch (via the atoms containing it)
E=set(S)
for t in cands:
    for a in L.var_atoms[t]:
        E |= set(L.atom2eq.get(a,{}))
E=sorted(E)
print(f"S={len(S)} failing ; constrained system over {len(E)} equations ; {len(cands)} candidate handles", flush=True)
def inner(vv, idx):
    a2=L.all_atom_values(vv)
    return [sum(c*a2[a] for a,c in L.eq_atoms[e][2].items()) for e in idx]
base=inner(v,E)
print("  nonzero inner sums at base:", sum(1 for x in base if x), flush=True)
cols=[];used=[]
t0=time.time()
for t in cands:
    old=v[t]
    v[t]=old+1; fw.forward(v); e1=inner(v,E)
    v[t]=old+2; fw.forward(v); e2=inner(v,E)
    v[t]=old; fw.forward(v)
    d1=[e1[i]-base[i] for i in range(len(E))]
    d2=[e2[i]-base[i] for i in range(len(E))]
    if all(d2[i]==2*d1[i] for i in range(len(E))) and any(d1):
        cols.append(d1); used.append(t)
print(f"  exact-linear handles over the FULL region: {len(used)} -> {used} ({time.time()-t0:.0f}s)", flush=True)
if used:
    M=[[cols[j][i] for j in range(len(used))] for i in range(len(E))]
    rhs=[-base[i] for i in range(len(E))]
    x=solve_int(M,rhs)
    print("  constrained integer solution:", "FOUND" if x else "NONE", flush=True)
    if x:
        for j,t in enumerate(used): v[t]+=x[j]
        fw.forward(v)
        f=L.failing_eqs(L.all_atom_values(v)); b=fw.bad_checks(v)
        print(f"  AFTER: failing={len(f)} score={L.NEQ-len(f)} bad_checks={len(b)} {b[:10]}")
        json.dump({('x_%d'%i):str(v[i]) for i in range(L.NVARS)}, open(os.path.join(HERE,'data','realise3_named.json'),'w'))
    else:
        for drop in range(1,10):
            hit=None
            for combo in itertools.combinations(range(len(S)),drop):
                keep=[i for i in range(len(E)) if E[i] not in {S[c] for c in combo}]
                if solve_int([M[i] for i in keep],[rhs[i] for i in keep]) is not None:
                    hit=combo; break
            if hit:
                print(f"  best: failing {drop} -> score {L.NEQ-drop}  (drop {[S[c] for c in hit]})")
                break
