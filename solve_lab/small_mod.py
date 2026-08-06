import json
from sympy import factorint
p=2**256-2**32-977
mods=[6672769,13523997,6122989,9994531,7376877,9367949,537773]
print("factorization of the small moduli:")
for m in mods:
    print(f"  {m} = {factorint(m)}")
# relation to p-1?
print(f"\np-1 factors (partial): ")
# gcd with p-1
from math import gcd
for m in mods:
    print(f"  gcd({m}, p-1) = {gcd(m,p-1)}")
# are they factors of any pin coef or related?
pinrec=json.load(open('pinrec.json'))
coefs=set(coef for i,sel,tgt,const,coef,handle in pinrec if coef>1)
print(f"\nsmall moduli that appear as pin coefs: {[m for m in mods if m in coefs]}")
# do the small moduli divide (const_i - const_j) for pin constants? (structure)
consts=[const%p for i,sel,tgt,const,coef,handle in pinrec]
for m in mods[:3]:
    cnt=sum(1 for c in consts if c%m==0)
    print(f"  pin consts divisible by {m}: {cnt}/{len(consts)}")
