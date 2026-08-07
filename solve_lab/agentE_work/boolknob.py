"""Treat the boolean selector movers as 0/1 DECISIONS (exact re-propagation, not derivatives)
   and ask whether any of them reaches a10187 or a20212 with a residual delta coprime to p."""
import sys, json, math, re, time, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import engine as E, fast, harness as H
P=115792089237316195423570985008687907853269984665640564039457584007908834671663
TARGET=[10187,20212]
CLUSTER=[7389,10187,20212,20215,28647]
s={int(k):int(v) for k,v in json.load(open(sys.argv[1] if len(sys.argv)>1 else 'triple8_seed.json')).items()}
v0=E.forward(s); bad0=E.badatoms(v0)
print("state bad:",sorted(bad0),flush=True)
def isb(f):
    for i in H.occ[f]:
        t=re.sub(r'x_%d\b'%f,'X',H.atoms[i])
        if t in ('X - X * X','X * X - X','X * (X - 1)','2 * X * (1 - X)'): return True
    return False
cand=set()
for a in TARGET: cand|=set(E.cone(a)[1])
cand=sorted(cand)
bools=[f for f in cand if isb(f)]
print(f"candidates in cones of {TARGET}: {len(cand)}; boolean: {len(bools)}",flush=True)
res={}; hits=[]
t0=time.time()
for f in bools:
    cur=v0[f]
    if cur not in (0,1):   # not currently a clean 0/1 -> still probe 0 and 1
        pass
    row={}
    for val in (0,1):
        if val==cur: 
            row[val]=bad0; continue
        b,_=fast.resid_delta(v0,bad0,{f:val})
        row[val]=b
    d={}
    for a in CLUSTER:
        d[a]=row[1].get(a,0)-row[0].get(a,0)
    coll=sorted((set(row[1])|set(row[0]))-set(bad0))
    res[f]=(d,coll,cur)
    for a in TARGET:
        if d[a] and math.gcd(abs(d[a]),P)==1:
            hits.append((f,a,d[a]))
print(f"probed {len(bools)} booleans in {time.time()-t0:.0f}s",flush=True)
print(f"*** boolean flips with delta COPRIME to p on {TARGET}: {len(hits)} ***",flush=True)
for f,a,d in hits[:40]:
    print(f"   x_{f} -> a{a}: delta bits={d.bit_length()} delta mod p = {d%P}",flush=True)
# summary of delta residues
import collections
cnt=collections.Counter()
for f,(d,coll,cur) in res.items():
    for a in TARGET:
        if d[a]==0: cnt[(a,'zero')]+=1
        elif d[a]%P==0: cnt[(a,'0 mod p')]+=1
        else: cnt[(a,'nonzero mod p')]+=1
print("delta census:",dict(cnt),flush=True)
pickle.dump(res,open('boolknob.pkl','wb'))
