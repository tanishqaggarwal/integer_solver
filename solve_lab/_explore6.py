import heal_harness as H, json
p=H.p
d=H.loadd('best/new_instance_partial_39013.json')
def setfree(dd):
    for v in range(H.NVARS): H.val[v]=dd.get(v,0)
setfree(d); H.forward()
F0=set(H.fails())
print(f"baseline 39013: {len(F0)} fail: {sorted(F0)}")

# Regime 1: set x_29322=0 and x_3558=0 exactly.
# x_29322 = x_14853 - x_12186 ; set x_14853 := x_12186
# x_3558 = x_24908 - x_16742 ; x_24908 is a gate, set x_16742 := forward value of x_24908
x24908 = H.val[24908]  # forward value (independent of the controls)
d2=dict(d)
d2[14853]=d[12186]          # x_29322 -> 0 exactly
d2[16742]=x24908            # x_3558 -> 0 exactly
setfree(d2); H.forward()
F=set(H.fails())
print(f"\nregime1 (x_29322=0,x_3558=0 exactly): {len(F)} fail")
print(f"  x_29322={H.val[29322]}  x_3558={H.val[3558]}")
print(f"  S=x_35389={H.val[35389]}  T=x_6671={H.val[6671]}")
print(f"  L1=x_11150={H.val[11150]}  L2=x_25739={H.val[25739]}  L3=x_37758={H.val[37758]}")
print(f"  fixed from F0: {sorted(F0-F)}")
print(f"  new fails: {sorted(F-F0)}")
print(f"  still failing (in both): {sorted(F0&F)}")
