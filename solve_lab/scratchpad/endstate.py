import sys, os
os.chdir('/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/scratchpad')
import heal_harness as H
import atomlib as A
p = H.p

vA = H.loadd('best_agentA_39022.json')
for v in H.freeinp:
    H.val[v] = vA.get(v, 0)
H.forward()

def zero_atom(ai):
    poly = A.ATOMS[ai]
    val = A.eval_atom(ai, H.val)
    for vl, c in poly:
        if len(vl) == 1 and vl[0] in H.freeinp:
            f = vl[0]
            if val % c == 0:
                H.val[f] -= val // c
                return f
    return None

# topological heal
for step in range(30):
    sn = [ai for ai in range(A.NATOM) if len(A.ATOM_VARS[ai]) <= 6 and A.eval_atom(ai, H.val) != 0]
    if not sn: break
    fixed = False
    for ai in sn:
        if zero_atom(ai) is not None: fixed = True
    H.forward()
    if not fixed: break

print(f"after topo heal: {len(H.fails())} fails; remaining small gaps:")
for ai in range(A.NATOM):
    if len(A.ATOM_VARS[ai]) <= 6 and A.eval_atom(ai, H.val) != 0:
        print(f"  atom {ai}: {A.ATOM_REPR[ai][:70]}  val%p={'0' if A.eval_atom(ai,H.val)%p==0 else 'NZ'}")

# crux residues
def r(x): return H.val[x] % p
print("\n=== CORE residues (mod p) ===")
print(f" x_14853 % p = {r(14853)}")
print(f" x_12186 % p = {r(12186)}")
print(f" x_14853 == x_12186 mod p ? {r(14853)==r(12186)}")
print(f" x_24908 % p = {r(24908)}")
print(f" x_16742 % p = {r(16742)}")
print(f" x_24908 == x_16742 mod p ? {r(24908)==r(16742)}")
print(f" x_29322 = x_14853 - x_12186 (%p) = {(H.val[14853]-H.val[12186])%p}")
print(f" x_3558  = x_24908 - x_16742 (%p) = {(H.val[24908]-H.val[16742])%p}")
print(f" x_1308 % p = {r(1308)}  (x_14853 target)")
print(f" x_31339 % p = {r(31339)}, x_6858 % p = {r(6858)} (x_31339 target)")

# core base gates S=x_35389, T=x_6671
print(f"\n S = x_35389 % p = {r(35389)}")
print(f" T = x_6671  % p = {r(6671)}")
