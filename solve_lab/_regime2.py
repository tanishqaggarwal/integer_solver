import heal_harness as H, json, pickle
p=H.p
d=H.loadd('best/new_instance_partial_39013.json')
def setfree(dd):
    for v in range(H.NVARS): H.val[v]=dd.get(v,0)
setfree(d); H.forward()
# baseline residues
a0=H.val[29322]%p          # x_29322
b0=H.val[3558]%p           # x_3558
x24453=H.val[24453]%p
x12186=H.val[12186]%p
x14853=H.val[14853]%p
x16742=H.val[16742]%p
x24908=H.val[24908]%p
print(f"a0={a0}\nb0={b0}")
# S=0: x_33469*a0^2 = b0^2 -> x_33469_target = b0^2 * inv(a0^2)
inv=lambda z: pow(z,p-2,p)
x33469_target=(b0*b0)%p * inv(a0*a0%p) %p
# x_33469 = x_22162 + x_12186 + x_14853 + x_24453
x22162_target=(x33469_target - x12186 - x14853 - x24453)%p
# x_1326 = x_12186 - x_22162 (after setting x_22162)
x1326_new=(x12186 - x22162_target)%p
# T=0: x_27713*a0 = b0*x_1326 -> x_27713_target = b0*x_1326/a0
x27713_target=(b0*x1326_new)%p * inv(a0) %p
# x_27713 = x_30213 + x_16742
x30213_target=(x27713_target - x16742)%p
print(f"x33469_target={x33469_target}")
print(f"x22162_target={x22162_target}")
print(f"x27713_target={x27713_target}")
print(f"x30213_target={x30213_target}")
# Set as integers >=0 congruent to target
d2=dict(d)
d2[22162]=x22162_target
d2[30213]=x30213_target
setfree(d2); H.forward()
print("\nAfter setting x_22162,x_30213:")
print(f"  x_29322%p={H.val[29322]%p} (a, should=a0)")
print(f"  x_3558%p={H.val[3558]%p} (b, should=b0)")
print(f"  x_33469%p={H.val[33469]%p} (target {x33469_target})")
print(f"  x_27713%p={H.val[27713]%p} (target {x27713_target})")
print(f"  S=x_35389%p={H.val[35389]%p}")
print(f"  T=x_6671%p={H.val[6671]%p}")
print(f"  L1%p={H.val[11150]%p}")
print(f"  L2%p={H.val[25739]%p}")
print(f"  L3%p={H.val[37758]%p}")
F=sorted(H.fails())
print(f"\n  fails: {len(F)}: {F}")
# save d2
json.dump({f"x_{k}":str(v) for k,v in d2.items()},open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/regime2.json','w'))
