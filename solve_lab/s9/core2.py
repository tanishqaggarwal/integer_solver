"""Can the second core's control differences u'=x_18123, w'=x_17576 be zeroed?"""
import pickle, sys
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
P=2**256-2**32-977
K1=33472904810391811973223207617762334363023286939839396241234196646906030803538671321618319
NV=38748
freeinp=[x for x in range(NV) if x not in definer]
v0=H.load_assignment('../best/new_instance_partial_39022.json')
ripple(v0,{2527:1,5096:K1,33612:0}); ripple(v0,{14853:v0[12186]})
ripple(v0,{7068:v0[2099]+7376877*v0[642],4432:v0[19964]+v0[28730]}); ripple(v0,{24548:v0[25442]})
print('state: u2=x_18123 =',v0[18123]%P)
print('       w2=x_17576 =',v0[17576]%P)
print('       x_30454=%d  x_10261=%d  x_16787=%d  x_25199=%d'%(v0[30454]%P,v0[10261]%P,v0[16787]%P,v0[25199]%P))
res={18123:[],17576:[]}
for f in freeinp:
    v=list(v0); ch,_=ripple(v,{f:v0[f]+1})
    for t in (18123,17576):
        if t in ch: res[t].append((f, (v[t]-v0[t])%P))
for t in (18123,17576):
    modp=[(f,d) for f,d in res[t] if d!=0]
    print(f'\nfree inputs moving x_{t}: {len(res[t])} total, {len(modp)} of them mod p')
    for f,d in modp[:12]: print(f'   x_{f}: d%p={d}')
