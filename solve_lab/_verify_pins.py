import heal_harness as H, json, pickle
p=H.p
d0=H.loadd('best/new_instance_partial_39013.json')
for v in range(H.NVARS): H.val[v]=d0.get(v,0)
H.forward()
print("Selectors (quadrant 1,1):")
for v in [7715,34554,15298,23597,19271,34606,5647,13913,25538,608,22978,38085,24530]:
    print(f"  x_{v} = {H.val[v]%p}   (==0:{H.val[v]==0})")
# CONST1, CONST2 from atoms
C=pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/atomcache.pkl','rb'))
atoms=C['atoms']
# atom1465 = x_24468 - CONST1 - x_32989 ; CONST1 is the literal
def getconst(ai):
    # find the [] term (constant)
    for vl,c in atoms[ai]['poly']:
        if len(vl)==0: return c
    return None
CONST1=-getconst(1465); CONST2=-getconst(602)
print(f"\nCONST1 (atom1465) mod p = {CONST1%p if CONST1 else None}")
print(f"x_22162 mod p        = {H.val[22162]%p}")
print(f"x_13682 mod p        = {H.val[13682]%p}")
print(f"x_22162 == CONST1 mod p? {H.val[22162]%p == (CONST1%p if CONST1 else -1)}")
print(f"x_13682 == CONST1 mod p? {H.val[13682]%p == (CONST1%p if CONST1 else -1)}")
# Now the key: are x_1308,x_23927,x_24908,x_19083 pinned? Search for pin-atoms (check atoms) on each
print("\n=== searching pin/check atoms on second-layer vars ===")
# a 'check' atom on var v: an atom containing v that is NOT the definer gate and has few free vars
gatedef=set(H.definer.keys())
for cv in [1308,23927,24908,19083]:
    print(f"\n-- x_{cv} (gate, def={H.gates[H.definer[cv]][1][:40] if cv in H.definer else 'FREE'}) --")
    cnt=0
    for ai,a in enumerate(atoms):
        vs=set()
        for vl,c in a['poly']: vs|=set(vl)
        if cv in vs and a['n_eq']>=8:  # real check atoms have higher n_eq
            cnt+=1
            if cnt<=6: print(f"    atom#{ai} n_eq={a['n_eq']}: {a['repr'][:75]}")
