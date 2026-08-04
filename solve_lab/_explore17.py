import heal_harness as H, json, pickle
p=H.p
d=H.loadd('best/new_instance_partial_39013.json')
for v in range(H.NVARS): H.val[v]=d.get(v,0)
H.forward()
print("=== wire multipliers in handles ===")
for v in [38100,11360,26064,15753,8972,10603]:
    print(f"  x_{v} = {H.val[v]}   (==p? {H.val[v]==p}) (==1? {H.val[v]==1})  free_anc={sorted(H.anc[v])}")
# what defines x_26064?
print("\nx_26064 definer:", H.gates[H.definer[26064]][1] if 26064 in H.definer else "FREE")
# trace x_26064 to its root
def trace(v,depth=0,maxd=8):
    ind="  "*depth
    if v in H.freeinp:
        print(f"{ind}x_{v}=FREE={H.val[v]}"); return
    gi=H.definer[v]; t,rhs,vids=H.gates[gi]
    print(f"{ind}x_{v} = {rhs}  (val={H.val[v]}, ==p:{H.val[v]==p})")
    if depth<maxd:
        for u in vids: trace(u,depth+1,maxd)
trace(26064,maxd=10)
print("\n=== condition x_24908 vs x_19083 mod p ===")
print(f"  x_24908 % p = {H.val[24908]%p}")
print(f"  x_19083 % p = {H.val[19083]%p}")
print(f"  diff mod p = {(H.val[24908]-H.val[19083])%p}")
# shared free ancestors
print(f"  x_24908 free_anc count={len(H.anc[24908])}, x_19083={len(H.anc[19083])}, shared={len(H.anc[24908]&H.anc[19083])}")
print(f"  x_24908-only anc: {sorted(H.anc[24908]-H.anc[19083])[:20]}")
