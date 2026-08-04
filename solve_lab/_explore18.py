import heal_harness as H, json, pickle
p=H.p
C=pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/atomcache.pkl','rb'))
atoms=C['atoms']; eq2atoms=C['eq2atoms']
d=H.loadd('best/new_instance_partial_39013.json')
for v in range(H.NVARS): H.val[v]=d.get(v,0)
H.forward()
# eq 2071 structure
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
print("=== eq 2071 ===")
print(lines[2071][:300])
print("\n atoms in eq2071:", eq2atoms.get(2071))
# M-handle definitions
for v in [4007,29804,35605,24453,30213,33469,27713,1326,22162]:
    df = H.gates[H.definer[v]][1] if v in H.definer else "FREE"
    print(f"x_{v} = {df}  free={v in H.freeinp} val%p={H.val[v]%p if v not in [4007,29804,35605] else H.val[v]}")
# does x_11150 need to be exactly 0? check the load values now (baseline, loads nonzero)
print("\nbaseline loads: L1=x_11150=",H.val[11150], " L2=x_25739=",H.val[25739]," L3=x_37758=",H.val[37758])
print("x_4007=",H.val[4007]," x_29804=",H.val[29804]," x_35605=",H.val[35605])
# verify L1 = 8646263*S+1073965*T exactly?
S=H.val[35389]; T=H.val[6671]
print("\n8646263*S+1073965*T =", 8646263*S+1073965*T, " vs L1=",H.val[11150], " equal:",8646263*S+1073965*T==H.val[11150])
print("10159099*S+6926539*T =", 10159099*S+6926539*T, " vs L2=",H.val[25739]," equal:",10159099*S+6926539*T==H.val[25739])
print("8272701*S+5921311*T =", 8272701*S+5921311*T, " vs L3=",H.val[37758]," equal:",8272701*S+5921311*T==H.val[37758])
