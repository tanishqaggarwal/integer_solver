import json, math
p=2**256-2**32-977
pinrec=json.load(open('pinrec.json'))
consts=sorted(set(const%p for i,sel,tgt,const,coef,handle in pinrec))
print(f"{len(consts)} distinct constants mod p")
g=0
for cc in consts: g=math.gcd(g,cc)
print(f"gcd of all: {g}")
small=[cc for cc in consts if cc<10**60 or cc>p-10**60]
print(f"constants <1e60 or >p-1e60: {len(small)}")
diffs=sorted(consts[i+1]-consts[i] for i in range(len(consts)-1))
print(f"smallest gap between sorted consts: {len(str(diffs[0]))} digits")
factms=[6672769,13523997,6122989,9994531,7376877,9367949,537773]
for cc in consts[:6]:
    fs=[m for m in factms if cc%m==0]
    print(f"  const ..{str(cc)[-8:]} divisible by: {fs}")
coefs=sorted(set(coef for i,sel,tgt,const,coef,handle in pinrec if coef>1))
print(f"distinct pin coefs ({len(coefs)}): {coefs[:25]}")
# check: are constants roots of small polys / related by the small coefs?
# ratio structure: c_i / c_j mod p small?
