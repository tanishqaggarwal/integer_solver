import heal_harness as H, json, pickle
p=H.p
d4=H.loadd('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/regime4.json')
for v in range(H.NVARS): H.val[v]=d4.get(v,0)
H.forward()
def trace(v,depth=0,maxd=3,seen=None):
    if seen is None: seen=set()
    ind="  "*depth
    if v in H.freeinp:
        print(f"{ind}x_{v}=FREE #anc=1"); return
    if v in seen or depth>maxd:
        print(f"{ind}x_{v}=... #anc={len(H.anc[v])}"); return
    seen.add(v)
    gi=H.definer[v]; t,rhs,vids=H.gates[gi]
    print(f"{ind}x_{v} = {rhs}   #anc={len(H.anc[v])} val==p:{H.val[v]==p}")
    for u in vids: trace(u,depth+1,maxd,seen)
print("###### x_13913 (C1 compensator candidate) ######")
trace(13913,maxd=3)
print("free anc:", sorted(H.anc[13913]))
print("\n###### x_25538 (C2 compensator candidate) ######")
trace(25538,maxd=3)
print("free anc:", sorted(H.anc[25538]))
# Do these share free ancestors with S/T cone or the controls?
Score=set()
for v in [33469,27713,29322,3558,1326,1308,19083,24908,23927]:
    Score|=H.anc[v]
print("\nx_13913 anc ∩ (S/T/pin cone):", sorted(H.anc[13913]&Score))
print("x_25538 anc ∩ (S/T/pin cone):", sorted(H.anc[25538]&Score))
print("x_13913 anc ∩ x_25538 anc:", sorted(H.anc[13913]&H.anc[25538]))
