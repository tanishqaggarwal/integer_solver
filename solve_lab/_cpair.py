import heal_harness as H, json, pickle
p=H.p
d4=H.loadd('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/regime4.json')
for v in range(H.NVARS): H.val[v]=d4.get(v,0)
H.forward()
def info(v):
    df=H.gates[H.definer[v]][1] if v in H.definer else "FREE"
    wire = (H.val[v]==p)
    print(f"  x_{v}: {'FREE' if v in H.freeinp else 'gate='+df[:40]}  val==p:{wire}  #anc={len(H.anc.get(v,[]))}")
print("C1 pair (x_24468):")
for v in [24468,32989,22399,11436,13682,34243]: info(v)
print("C2 pair (x_18956):")
for v in [18956,14257,23917,7497,37892,32237]: info(v)
# Do x_22399, x_23917 == p (wire)? then x_32989=p*x_11436, x_14257=p*x_7497
print("\nx_22399==p:",H.val[22399]==p," x_23917==p:",H.val[23917]==p)
print("x_11436 free:",11436 in H.freeinp," x_7497 free:",7497 in H.freeinp)
print("x_34243 free:",34243 in H.freeinp," x_32237 free:",32237 in H.freeinp)
# does x_13682 actually depend on x_22162 numerically? and x_37892 on x_30213?
d0=H.loadd('best/new_instance_partial_39013.json')
for v in range(H.NVARS): H.val[v]=d0.get(v,0)
H.forward()
b13682=H.val[13682]; b37892=H.val[37892]
for v in range(H.NVARS): H.val[v]=d0.get(v,0)
H.val[22162]=d0.get(22162,0)+p  # change by p
H.forward()
print(f"\nchange x_22162 by p: x_13682 delta = {(H.val[13682]-b13682)}")
for v in range(H.NVARS): H.val[v]=d0.get(v,0)
H.val[30213]=d0.get(30213,0)+p
H.forward()
print(f"change x_30213 by p: x_37892 delta = {(H.val[37892]-b37892)}")
