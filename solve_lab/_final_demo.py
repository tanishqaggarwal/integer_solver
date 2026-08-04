import heal_harness as H, pickle
p=H.p
C=pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/atomcache.pkl','rb'))
atoms=C['atoms']
def av(ai):
    a=atoms[ai];s=0
    for vl,c in a['poly']:
        t=c
        for v in vl:t*=H.val[v]
        s+=t
    return s
d0=H.loadd('best/new_instance_partial_39013.json')
def setb():
    for v in range(H.NVARS): H.val[v]=d0.get(v,0)
setb();H.forward()
F0=set(H.fails()); S0=H.val[35389]%p;T0=H.val[6671]%p
inv=lambda z:pow(z,p-2,p)
print("=== 1-DOF demonstration (quadrant 1,1) ===")
print(f"baseline: S={S0!=0 and 'nonzero' or 'ZERO'}, T={T0!=0 and 'nonzero' or 'ZERO'}, {len(F0)} fails")
# clean T-knob x_30213: solve T=0 -> x_27713 = b0*x1326/a0 ; x_27713=x_30213+x_16742
a0=H.val[29322]%p;b0=H.val[3558]%p;x1326=H.val[1326]%p;x16742=H.val[16742]%p
x27713_t=(b0*x1326)%p*inv(a0)%p
x30213_t=(x27713_t-x16742)%p
setb();H.val[30213]=x30213_t;H.forward()
F=set(H.fails())
print(f"\napply clean T-knob (x_30213): T={H.val[6671]%p!=0 and 'nonzero' or 'ZERO'}, S={H.val[35389]%p==S0 and 'UNCHANGED(conserved)' or 'moved'}")
print(f"  wiring breaks beyond core: {sorted(F-F0)}  (empty=clean); total fails={len(F)}")
# Now try to move S: only x_22162 does it directly -> breaks C1 (atom29373)
x33469_t=(b0*b0)%p*inv(a0*a0%p)%p
x22162_t=(x33469_t-H.val[12186]%p-H.val[14853]%p-H.val[24453]%p)%p
setb();H.val[22162]=x22162_t;H.forward()
print(f"\napply S-knob (x_22162): S={'ZERO' if H.val[35389]%p==0 else 'moved'}")
print(f"  but breaks guard 9648 (C1 pin atom29373)={av(29373)%p!=0}; C1 is p-granular(x_9254=p*x_33787,x_34243=p*x_14393) -> UNHEALABLE")
print(f"  fails now={len(H.fails())}")
print("\n=> S can ONLY be moved by breaking the C1 message-load pin (x_22162==CONST1 mod p), which is p-granular and cannot be reheal ed. S is conserved. DOF=1.")
