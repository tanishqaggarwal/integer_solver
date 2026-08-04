#!/usr/bin/env python3
import heal_harness as H
p=H.p
F16=[697,1985,5225,10815,16048,17784,17801,22402,23667,24721,27124,28737,29638,29959,35935,37431]
# free-ancestor sets of the 16
anc16=set()
for i in F16:
    for w in H.eqvars[i]:
        if w in H.freeinp: anc16.add(w)
        anc16|=H.anc.get(w,set())
# core knobs and whether they're in anc16
knobs={'x_14853':14853,'x_12186':12186,'x_16742':16742,'x_31339':31339}
for lbl,k in knobs.items():
    print(f"{lbl}: feeds the 16? {k in anc16}")
# Also, does the CORE (x_29322,x_3558) depend on knobs that avoid the 16?
# x_29322 = x_14853-x_12186; x_3558 = x_24908-x_16742. x_24908 free-ancestors:
anc24908=H.anc.get(24908,set())
avoid=[w for w in anc24908 if w not in anc16]
print(f"\nx_24908 free-ancestors NOT feeding the 16: {len(avoid)}/{len(anc24908)}")
# which core-control free inputs are OUTSIDE anc16 (can move without breaking the 16)?
core_ctrl = {14853,12186,16742} | anc24908
safe=[w for w in core_ctrl if w not in anc16 and w in H.freeinp]
print(f"core-control free inputs that DON'T feed the 16: {sorted(safe)[:20]} (count {len(safe)})")
