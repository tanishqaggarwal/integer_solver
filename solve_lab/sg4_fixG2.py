import heal_harness as H
p=H.p
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
V0=H.val[:]
baseF=set(H.fails())
print(f"baseline: {len(baseF)} fails: {sorted(baseF)}")
# Fix G2 via x_12553 (shifts x_19964 mod p) + x_9413 quotient
R2=(H.val[19964]-H.val[4432])%p
H.val[12553]=V0[12553]-R2   # now x_19964 = x_12553 ≡ x_4432 mod p
H.forward()
d=H.val[4432]-H.val[19964]
print("after x_12553-=R2: (x_4432-x_19964) =",d," divisible by p?",d%p==0)
if d%p==0:
    H.val[9413]=d//p
    H.forward()
def ev642(): 
    return sum(c*(H.val[m[0]]*H.val[m[1]] if len(m)==2 else (H.val[m[0]] if m else 1)) for m,c in [((),0)])
# check G2
g2=H.val[4432]-H.val[19964]-H.val[28730]
print("G2 =",g2)
F=set(H.fails())
print(f"after fixing G2: {len(F)} fails: {sorted(F)}")
print("  fixed:",sorted(baseF-F))
print("  broke:",sorted(F-baseF))
