import heal_harness as H, json, random
p=H.p
d=H.loadd('best/new_instance_partial_39013.json')
for v in range(H.NVARS): H.val[v]=d.get(v,0)
H.forward()
F0=set(H.fails())
print(f"baseline 39013: {len(F0)} fail")
print("control frees current values:")
for v in [14853,12186,16742]:
    print(f"  x_{v} = {H.val[v]}  (is free: {v in H.freeinp})")
print(f"  x_24908(gate) = {H.val[24908]}")
print(f"  x_29322 = x_14853-x_12186 = {H.val[29322]}  mod p = {H.val[29322]%p}")
print(f"  x_3558 = x_24908-x_16742 = {H.val[3558]}  mod p = {H.val[3558]%p}")

# Perturb x_14853 by a random amount, re-forward, check fails
import copy
base=[H.val[v] for v in range(H.NVARS)]
for delta in [1, 12345, p]:
    for v in range(H.NVARS): H.val[v]=d.get(v,0)
    H.val[14853]=d.get(14853,0)+delta
    H.forward()
    F=set(H.fails())
    print(f"perturb x_14853 += {delta}: {len(F)} fail; new fails not in F0: {sorted(F-F0)[:20]} ; fixed: {sorted(F0-F)}")
