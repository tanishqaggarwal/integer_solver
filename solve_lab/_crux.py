import heal_harness as H
p=H.p
d=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward(); V=H.val
# dependency: are x_2099, x_19964 free? do they depend on the gap knobs?
knobs=[7068,4432,17325,9413,2099,19964]
for v in [2099,19964,7068,4432,17325,9413]:
    isfree = v in H.freeinp
    anc = H.anc.get(v,set())
    dep = sorted((anc & set(knobs)) - {v})
    print(f"x_{v}: {'FREE' if isfree else 'GATE'}  depends on knobs {dep}")
# G1 divisibility: is (x_7068 - x_2099) divisible by 7376877*p?
g1const = V[7068]-V[2099]
M1=7376877*p
print(f"\nG1: need 7376877*p*x_17325 = x_7068 - x_2099 = {g1const}")
print(f"   7376877*p = {M1}")
print(f"   (x_7068-x_2099) mod (7376877*p) = {g1const % M1}")
print(f"   divisible? {g1const % M1 == 0}  quotient={g1const//M1 if g1const%M1==0 else 'N/A'}")
# also mod p and mod 7376877 separately
print(f"   (x_7068-x_2099) mod p = {g1const % p}")
print(f"   (x_7068-x_2099) mod 7376877 = {g1const % 7376877}")
# G2 divisibility: is (x_4432 - x_19964) divisible by p?
g2const=V[4432]-V[19964]
print(f"\nG2: need p*x_9413 = x_4432 - x_19964 = {g2const}")
print(f"   (x_4432-x_19964) mod p = {g2const % p}")
print(f"   divisible by p? {g2const % p == 0}  quotient={g2const//p if g2const%p==0 else 'N/A'}")
