import heal_harness as H
p=H.p
# print gate definitions for key gates
def gdef(v):
    if v in H.freeinp: return f"x_{v}=FREE"
    gi=H.definer[v]; t,rhs,vids=H.gates[gi]
    return f"x_{v} = {rhs}"
for v in [29322,3558,24908,16742,14853,12186,33469,29356,27713,1326,35389,6671,32680,11602,27762,11150,25739,37758]:
    print(gdef(v))
print("---- free anc of S=x_35389 ----")
print(sorted(H.anc[35389]))
print("---- free anc of T=x_6671 ----")
print(sorted(H.anc[6671]))
print("---- free anc of x_24908 ----")
print(sorted(H.anc[24908]))
print("---- union all core loads free anc ----")
U=set()
for v in [11150,25739,37758,35389,6671,3558,29322]:
    U|=H.anc[v]
print("size",len(U))
print(sorted(U))
