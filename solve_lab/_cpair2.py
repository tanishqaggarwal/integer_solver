import heal_harness as H, json, pickle
p=H.p
d4=H.loadd('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/regime4.json')
for v in range(H.NVARS): H.val[v]=d4.get(v,0)
H.forward()
def info(v):
    df=H.gates[H.definer[v]][1] if v in H.definer else "FREE"
    print(f"  x_{v}: {'FREE' if v in H.freeinp else 'gate='+df[:45]}  val==p:{H.val[v]==p}  val==1:{H.val[v]==1}  #anc={len(H.anc.get(v,[]))} anc={sorted(H.anc.get(v,[]))[:6]}")
print("x_34243 = x_16153 * x_14393:")
for v in [34243,16153,14393]: info(v)
print("x_32237 = x_21023 * x_22820:")
for v in [32237,21023,22820]: info(v)
# what free inputs control x_34243 and x_32237 residues?
print("\nx_34243 free anc:", sorted(H.anc[34243]))
print("x_32237 free anc:", sorted(H.anc[32237]))
