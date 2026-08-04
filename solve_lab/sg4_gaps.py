import heal_harness as H
p=H.p
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
V=H.val
print("=== G1/G2 slack structure at agentA baseline ===")
print("x_642   =", V[642], " = x_28599*x_17325 ; x_28599=",V[28599]%p==0 and 'p-mult' or V[28599], " x_17325=",V[17325])
print("x_28599 mod p =", V[28599]%p, " (0 => wire=p)")
print("x_28730 =", V[28730], " = x_17499*x_9413 ; x_17499 mod p=",V[17499]%p, " x_9413=",V[9413])
print()
G1 = 7376877*V[642] + V[2099] - V[7068]
G2 = V[4432] - V[19964] - V[28730]
print("G1 = 7376877*x_642 + x_2099 - x_7068 =", G1)
print("   x_642=",V[642]," x_2099=",V[2099]," x_7068=",V[7068])
print("G2 = x_4432 - x_19964 - x_28730 =", G2)
print("   x_4432=",V[4432]," x_19964=",V[19964]," x_28730=",V[28730])
print()
print("=== What G1=0 needs ===")
print("7376877*x_642 = x_7068 - x_2099 =", V[7068]-V[2099])
print("  (x_7068-x_2099) / 7376877 =", (V[7068]-V[2099])/7376877, " exact?", (V[7068]-V[2099])%7376877==0)
print("=== What G2=0 needs ===")
print("x_28730 = x_4432 - x_19964 =", V[4432]-V[19964])
print("  need x_9413 = that / x_17499 ; x_17499=",V[17499])
d=V[4432]-V[19964]
print("  (x_4432-x_19964) mod p =", d%p, " => needs x_28730 sub-p, but x_17499 mult of p")
