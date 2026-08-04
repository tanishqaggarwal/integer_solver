import heal_harness as H
p=H.p
print("NVARS",H.NVARS,"n gates",len(H.gates),"n freeinp",len(H.freeinp),"n eqs",len(H.eqcode))
# key vars
core_targets=[35389,6671,11150,25739,37758,3558,29322,14853,12186,24908,16742,33469,29356,27713,1326,32680,11602,27762]
for v in core_targets:
    a=H.anc.get(v,set())
    freea=a  # anc are free ancestors already
    print(f"x_{v}: #free_anc={len(freea)} {'FREE' if v in H.freeinp else 'gate'}")
