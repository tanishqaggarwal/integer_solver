"""Close G1/G2 (leaves=selected constants), then scan ripple-side selector combos.
Question: can the 16 ripple fails drop below 11?"""
import _bitlab as L, heal_harness as H, itertools, json, time
p=H.p
def apply_close(bits):
    """apply bit pattern, then close G1/G2 by setting leaves=x_2099,x_19964, forward again."""
    F=L.apply_pattern(bits,twopass=False)  # sets pins, forwards
    V=H.val
    # close G1: x_7068 = x_2099 + 7376877*x_642 ; x_642 currently 0 -> x_7068=x_2099
    V[7068]=V[2099]+7376877*V[642]
    # close G2: x_4432 = x_19964 + x_28730 ; x_28730 currently 0
    V[4432]=V[19964]+V[28730]
    H.forward()
    return H.fails()

base=dict(L.AGENTA_BITS)
# first: agentA + close G1/G2
F=apply_close(base)
print('agentA + close G1/G2:',len(F),'fails')
# ripple-side selectors
RSEL=[5910,11368,13195,17406,18022,22562,23751,28005]
res=[]
t=time.time()
for combo in itertools.product([0,1],repeat=8):
    b=dict(base)
    for s,v in zip(RSEL,combo): b[s]=v
    F=apply_close(b)
    res.append((len(F),combo))
res.sort(key=lambda x:x[0])
print('scanned 256 close-G1G2 combos in %.1fs'%(time.time()-t))
print('best 12:')
for n,combo in res[:12]:
    on=[RSEL[i] for i in range(8) if combo[i]]
    print(f'  {n} fails  ripple-ON={on}')
from collections import Counter
print('dist:',sorted(Counter(n for n,_ in res).items())[:10])
