import sys, os
os.chdir('/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
import heal_harness as H
p = H.p

vA = H.loadd('best_agentA_39022.json')
for v in H.freeinp:
    H.val[v] = vA.get(v, 0)
H.forward()
F0 = H.fails()
print(f"baseline fails: {len(F0)}: {F0}")

# read current computed x_2099, x_19964
x2099 = H.val[2099]
x19964 = H.val[19964]
x642 = H.val[642]
x28730 = H.val[28730]
print(f"x_2099={x2099%p} x_19964(%p)={x19964%p} x_642={x642} x_28730={x28730}")
print(f"G1 = x_7068 - x_2099 - 7376877*x_642 = {H.val[7068]-x2099-7376877*x642}")
print(f"G2 = x_4432 - x_19964 - x_28730 = {H.val[4432]-x19964-x28730}")

# THE OBVIOUS MOVE: set free leaves to satisfy G1, G2 (with x_642=x_28730=0 currently)
H.val[7068] = x2099 + 7376877 * x642
H.val[4432] = x19964 + x28730
H.forward()
F1 = H.fails()
print(f"\nafter x_7068:=x_2099+..., x_4432:=x_19964+...: {len(F1)} fails: {F1}")
newfails = set(F1) - set(F0)
fixed = set(F0) - set(F1)
print(f"  newly broken: {sorted(newfails)}")
print(f"  fixed: {sorted(fixed)}")
