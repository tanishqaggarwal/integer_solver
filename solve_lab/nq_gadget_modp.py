import heal_harness as H
p=H.p
d=H.loadd('sy_regime11_39018.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
H.forward()
V=H.val
R12553=42775533402728869434716629464193396056515231264222641773817154079369026410240838606908039
print(f"x_12553%p       = {V[12553]%p}")
print(f"R12553%p        = {R12553%p}")
print(f"match: {V[12553]%p==R12553%p}")
# gadget mod p analysis
x7181=V[7181]%p; x17925=V[17925]%p; x14865=V[14865]%p
print(f"\nx_7181%p = {x7181}")
# is x_7181 a QR mod p?
ls=pow(x7181,(p-1)//2,p)
print(f"Legendre(x_7181/p) = {ls}  ({'QR' if ls==1 else 'NON-residue' if ls==p-1 else 'zero'})")
# current gadget residuals
print(f"\ncurrent x_27177%p = {V[27177]%p}")
print(f"current x_4306%p  = {V[4306]%p}")
print(f"current x_31731%p = {V[31731]%p}")
# required x_27019 residue for x_27177=0: x_27019^2 = x_7181*x_17925^2
if ls==1:
    import sympy
    r=sympy.sqrt_mod(x7181,p)  # sqrt of x_7181
    x27019_needed=(r*x17925)%p
    x12553_needed1=(x14865+x27019_needed)%p
    x12553_needed2=(x14865-x27019_needed)%p
    print(f"\nsqrt(x_7181)%p = {r}")
    print(f"x_27019 needed (±) = {x27019_needed} or {(-x27019_needed)%p}")
    print(f"x_12553 needed = {x12553_needed1}")
    print(f"           or   = {x12553_needed2}")
    print(f"R12553%p       = {R12553%p}")
    print(f"match1={x12553_needed1==R12553%p} match2={x12553_needed2==R12553%p}")
