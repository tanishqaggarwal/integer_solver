import sys, os
os.chdir('/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/scratchpad')
import heal_harness as H
import atomlib as A
p = H.p

msg = [6418, 9118, 31861, 8731, 12553, 14865, 2081, 4287]
# descendants: gate outputs whose anc contains the msg input
for m in msg:
    desc = [t for t in H.order if m in H.anc[t]]
    natoms = len(A.VAR_ATOMS[m])
    role = 'FREE' if m in H.freeinp else 'GATE'
    print(f"x_{m}: {role}, in {natoms} atoms, feeds {len(desc)} gate-outputs")

# check: does x_2099 depend LINEARLY on x_6418,x_9118,x_31861? finite-diff test
vA = H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
base2099 = H.val[2099]; base19964 = H.val[19964]
print(f"\nbase x_2099%p={base2099%p}")
print(f"G1 gap = (x_7068 - x_2099)%p = {(H.val[7068]-H.val[2099])%p}  (need 0 for G1)")
print(f"G2 gap = (x_4432 - x_19964)%p = {(H.val[4432]-H.val[19964])%p}")

# perturb each msg input by +1 and see effect on x_2099, x_19964, and #fails
import copy
F0 = set(H.fails())
print(f"\nbaseline fails {len(F0)}")
for m in [6418, 9118, 31861, 8731, 12553, 14865]:
    old = H.val[m]
    H.val[m] = old + 1
    H.forward()
    d2099 = (H.val[2099]-base2099)%p
    d19964 = (H.val[19964]-base19964)%p
    F = set(H.fails())
    newbroken = F - F0
    print(f"  x_{m}+=1: d(x_2099)={d2099 if d2099< p//2 else d2099-p}, d(x_19964)={d19964 if d19964<p//2 else d19964-p}, fails={len(F)} newbroken={sorted(newbroken)[:8]}")
    H.val[m] = old
    H.forward()
