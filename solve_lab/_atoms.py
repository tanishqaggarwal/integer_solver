import heal_harness as H
p=H.p
d=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward(); V=H.val
def show(name,expr):
    print(f"  {name} = {expr}   (bits {abs(expr).bit_length()}, mod p={expr%p}, /p={'p*'+str(expr//p) if expr%p==0 else 'NOT mult of p'})")
print("Wire checks (should be 0 if wire correct):")
show("A1: x_642 - x_28599*x_17325", V[642]-V[28599]*V[17325])
show("A3: x_28730 - x_17499*x_9413", V[28730]-V[17499]*V[9413])
show("A5: x_23754 - x_26874*x_6947", V[23754]-V[26874]*V[6947])
print("wire values: x_28599=%d (=p? %s), x_17499=%d (=p? %s), x_26874=%d (=p? %s)"%(V[28599],V[28599]==p,V[17499],V[17499]==p,V[26874],V[26874]==p))
print("\nGap atoms:")
show("G1: 7376877*x_642 + x_2099 - x_7068", 7376877*V[642]+V[2099]-V[7068])
show("G2: x_4432 - x_19964 - x_28730", V[4432]-V[19964]-V[28730])
print("  x_642=%d (bits %d), x_2099=%d, x_7068=%d"%(V[642],V[642].bit_length(),V[2099],V[7068]))
print("  x_4432=%d, x_19964=%d, x_28730=%d (bits %d)"%(V[4432],V[19964],V[28730],V[28730].bit_length()))
print("\ngadget atoms (eq12270):")
show("x_29854 - x_22665*x_1329", V[29854]-V[22665]*V[1329])
show("5113045*x_7075*x_9118 - x_29854", 5113045*V[7075]*V[9118]-V[29854])
show("x_31864 - x_28961*x_10903", V[31864]-V[28961]*V[10903])
show("x_7075*x_8731 + x_31864", V[7075]*V[8731]+V[31864])
print("  x_7075=%d (=p? %s), x_22665=%d(=p?%s), x_28961=%d(=p?%s)"%(V[7075],V[7075]==p,V[22665],V[22665]==p,V[28961],V[28961]==p))
print("  x_1329=%d x_9118=%d x_8731=%d x_10903=%d x_29854=%d x_31864=%d"%(V[1329],V[9118],V[8731],V[10903],V[29854],V[31864]))
