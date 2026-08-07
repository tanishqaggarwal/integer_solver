"""Does flipping a selector change the a10187 residue multiset?  (The counting argument in
   LOG 15.1 assumes it is fixed.)"""
import sys, json, math, pickle, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import engine as E, fast
P=115792089237316195423570985008687907853269984665640564039457584007908834671663
res=pickle.load(open('boolknob.pkl','rb'))
base={int(k):int(v) for k,v in json.load(open('triple8_seed.json')).items()}
ALL=sorted(res)                       # the 256 boolean movers probed before
def census(seed, rows=(10187,20212)):
    v0=E.forward(seed); bad0=E.badatoms(v0)
    out={}
    for a in rows: out[a]=collections.Counter()
    for f in ALL:
        if seed.get(f,0)==1:  continue          # already on
        b1,_=fast.resid_delta(v0,bad0,{f:1})
        for a in rows:
            d=b1.get(a,0)-bad0.get(a,0)
            out[a][ (d%P) if d else 'zero' ]+=1
    return out, bad0
print("=== baseline (selectors 1530,1603 on) ===",flush=True)
c0,bad0=census(base)
for a in (10187,20212):
    print(f" a{a}: {len(c0[a])} classes; " + ", ".join(f"{str(k)[:22]}..x{v}" for k,v in c0[a].most_common(5)),flush=True)
R1=bad0.get(20215,0); R2=bad0.get(28647,0)
print(" R1 mod p =",R1%P," R2 mod p =",R2%P,flush=True)
for probe in (490, 2081, 4287, 5910, 12054):
    s=dict(base); s[probe]=1
    c,bd=census(s)
    print(f"=== after flipping x_{probe} ===  bad={sorted(bd)[:8]}",flush=True)
    for a in (10187,20212):
        same = (c[a]==c0[a])
        print(f" a{a}: {len(c[a])} classes, identical to baseline? {same}; "
              + ", ".join(f"{str(k)[:22]}..x{v}" for k,v in c[a].most_common(4)),flush=True)
    print(f"  R1 mod p = {bd.get(20215,0)%P}   R2 mod p = {bd.get(28647,0)%P}",flush=True)
