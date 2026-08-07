import sys, collections, pickle, time, json
sys.path.insert(0,'.')
import common as C
import harness as H, engine as E, fast
P=C.P
struct=pickle.load(open('struct_supp.pkl','rb'))
seed=dict(C.BASE); v0=E.forward(seed); bad0=E.badatoms(v0)
FOOT=collections.defaultdict(set)
for e,(issq,outer,terms) in enumerate(H.eqt):
    for c,a in terms:
        if a>=0: FOOT[a].add(e)
NF={a:len(s) for a,s in FOOT.items()}
pure=[f for f in struct if len(struct[f])==1]
print("pure-structural handles:",len(pure))
res={}; t0=time.time()
for f in pure:
    a=next(iter(struct[f])); o=v0[f]
    b1,_=fast.resid_delta(v0,bad0,{f:o+1}); b2,_=fast.resid_delta(v0,bad0,{f:o+2})
    d1=b1.get(a,0)-bad0.get(a,0); d2=b2.get(a,0)-bad0.get(a,0)
    res[f]=(a,d1,d2,d2==2*d1)
print("measured in %.0fs"%(time.time()-t0))
byatom=collections.defaultdict(list)
for f,(a,d1,d2,aff) in res.items(): byatom[a].append((f,d1,aff))
print("distinct atoms with a pure handle:",len(byatom))
# modulus per atom = gcd of affine steps
import math
mods={}
for a,lst in byatom.items():
    g=0
    for f,d1,aff in lst:
        if aff and d1: g=math.gcd(g,abs(d1))
    mods[a]=g
pk=collections.Counter()
for a,g in mods.items():
    if g==0: pk['zero-step']+=1
    elif g==P: pk['exactly p']+=1
    elif g%P==0: pk['multiple of p']+=1
    else: pk['other']+=1
print("handle modulus classes:",dict(pk))
oth=[(a,mods[a]) for a in mods if mods[a] and mods[a]%P]
print("non-p moduli (first 15):",[(a,mods[a].bit_length(),math.gcd(mods[a],P)) for a,_ in oth[:15]])
pickle.dump({'res':res,'byatom':dict(byatom),'mods':mods,'NF':NF},open('handles.pkl','wb'))
print("\nkey atoms:")
for a in (20215,28647,20212,7389,10187,747,30787,26958,40306,20649,20652,32148,28033,28035,28037):
    print(f"  a{a}: nf={NF.get(a)} handles={len(byatom.get(a,[]))} mod={'p' if mods.get(a)==P else (str(mods[a]//P)+'*p' if mods.get(a) and mods[a]%P==0 else mods.get(a))}")
